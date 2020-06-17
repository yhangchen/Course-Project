
                     # ------------- COLUMN GENERATION ALGORITHM FOR THE MEDIUM INSTANCE ------------- 

from __future__ import print_function

import sys
import time

import cplex
from cplex.exceptions import CplexSolverError

import numpy as np
np.set_printoptions(threshold=np.nan)

import openpyxl

from medium_data import getNodes, getArcs, getCommodities, getCosts # python script to get the data


def dijkstra(origin):
    """ Return the shortest paths from origin
        If mostUsedNode != -1, delete it from the graph, then compute the shortest paths from origin
        (The node will be considered as deleted for the next iterations too)
    """
#   -------- DIJKSTRA'S ALGORITHM --------

    if(mostUsedNode != -1):
        for i in range(nStations):
            cost[i, mostUsedNode, 0] = 0 # Disconnect the mostUsedNode from all others (cumulative)
                
    labels = np.full((nStations), np.inf) # Contains the cost of the shortest paths from origin with init value = infinity
    labels[origin] = 0
    precedent_nodes = np.full((nStations), -1, dtype=int) # Contains the previous node of each node in the shortest path
    toTreat = np.array([origin]) # Node that need to be processed
    
    while True:
        currentNode = toTreat[np.argmin(np.take(labels, toTreat))] # Pick the node with lowest shortest path from the set toTreat
        for i in range(nStations):
            if(cost[currentNode, i, 0] > 0): #Among all neighbors of currentNodes
                if(labels[i] > labels[currentNode] + cost[currentNode, i, 0]): # If dual constraint violated
                    
                    labels[i] = labels[currentNode] + cost[currentNode, i, 0] # Lower it to satisfy the feasibility and slackness property
                    precedent_nodes[i] = currentNode # Get the node_id
                    toTreat = np.append(toTreat, i) # As we found a better path to go to i, need check on the neighbors of i

        toTreat = np.delete(toTreat,np.where(toTreat==currentNode)[0]) # Delete the node we just processed
        if(toTreat.size == 0): # If there is no node left to treat, we are sure we obtained the best solution.
            break
        
    return labels, precedent_nodes

def getInitSet():
    """ Return a feasible set of inital variables.
        If mostUsedNode != -1, we look for the shortest paths without it
        Otherwise, just run Dijsktra's algorithm on all the graph
    """
    global shortest_paths, previous_nodes, pathsK_nodes

#   -------- GET INITIAL SET OF PROMISSING VARIABLES --------
    
    distance_done = set() # Nodes for which dijkstra's has already been run
    for i in range(nCommodities): #Don't need the last node, if you have already find all the shortest path from every node except the last, you already computed the info for the last node
        print(i)
        start = commodities[i][0]
        dest = commodities[i][1]
        if(mostUsedNode != -1):
            #Check if mostUsedNode is in the global shortest path and isn't the destination of the commodity (Otherwise, no need to recompute shorest path)
            if(sum([pathsK_nodes[i][0][index] for index in usedNode]) > 0 and commodities[i][1] not in usedNode):

                if(start not in distance_done): # If we haven't computed the shortest paths from the commoditiy origin, do it
                    shortest_paths[start], previous_nodes[start] = dijkstra(start)
            
                if(previous_nodes[start, dest] != -1): # Check if the mostUsedNode deletion hasn't disconnected the destination from the origin
                    
                    newPath_nodes = [0] * nStations
                    newPath_nodes[dest] = 1 * commodities[i][2]
                    # Create the new path going backward
                    while(dest != start):
                        newPath_nodes[previous_nodes[start, dest]] = 1 * commodities[i][2]
                        dest = previous_nodes[start, dest]
        
                    distance_done.add(start)

                    newPath_nodes[start] = 0 #Exclude the starting node
                    newPath_nodes.append( shortest_paths[start, commodities[i][1]] )
                    if(newPath_nodes not in pathsK_nodes[i]): # Check if the path is already taken into account
                        pathsK_nodes[i].append(newPath_nodes)
                        
        #Else : Just run Dijkstra's on whole graph
        else:
            if(start not in distance_done): # If we haven't computed the shortest path from the commoditiy origin yet, do it
                shortest_paths[start], previous_nodes[start] = dijkstra(start)

            newPath_nodes = [0] * nStations
            newPath_nodes[dest] = 1 * commodities[i][2]
                
            while(dest != start):
                newPath_nodes[previous_nodes[start, dest]] = 1 * commodities[i][2]
                dest = previous_nodes[start, dest]
        
            distance_done.add(start)
        
            newPath_nodes[start] = 0
            newPath_nodes.append( shortest_paths[start, commodities[i][1]] )
            pathsK_nodes[i].append(newPath_nodes)

    return pathsK_nodes

def getMostUsedNode(): 
    """ Return the node liable to transport the biggest quantity minus its capacity
        (Probabilities a path is taken are equiprobable)
    """
    count_node = np.zeros(nStations) # Contains the exceeding quantity for each node
    for i in range(nStations):
        count_node[i] = count_node[i] - capacity_node[node[i][2]-1]
        for j in range(nCommodities):
            total = len(pathsK_nodes[j])
            for k in range(len(pathsK_nodes[j])):
                count_node[i] += pathsK_nodes[j][k][i]/total # (Quantity of the commodity)*(Proba the path is choosen)
                
    return np.argmax(count_node)

def dual_of_Restrited(pathsK_nodes):
    """ Goal : Obtain the dual variables of the restricted master problem
    """
#   ------------- DUAL OF THE RESTRICTED MASTER PROBLEM -------------
    
    model = cplex.Cplex() # Initialize the model
    model.objective.set_sense(model.objective.sense.maximize) ## Set the objective function to maximization

    #Add variables
    # Correspond to the flow conservation constraints
    for i in range(nCommodities):
        model.variables.add(obj = [1], lb = [-cplex.infinity], ub = [cplex.infinity],#obj contains the coefficients of the decision variables in the objective function
                            types = [model.variables.type.continuous],
                            names = ['Y( K_' + str(i) + ')'])
        
    # Correspond to the node capacity constraints 
    for i in range(nStations):
        model.variables.add(obj = [ capacity_node[node[i][2]-1] ], lb = [-cplex.infinity], ub = [0],
                            types = [model.variables.type.continuous],
                            names = ['Y( ' + str(node[i][0])  + ')'])
        
    #Add constraints
    for i in range(nCommodities):
        for j in range(len(pathsK_nodes[i])):
            ind = []# Put the indices of the non-zero variables of the current constraint
            val = []# Put their coefficients here
            row = []
            ind.append(i)
            val.append(1)
            for k in range(nStations):
                if(pathsK_nodes[i][j][k] > 0): #Check which node is contained in the current paths (do not include starting node)
                    ind.append(nCommodities+k)
                    val.append(commodities[i][2])
            row.append([ind, val])
            model.linear_constraints.add(lin_expr = row,
                                     senses = "L", #Less or equal than constraint
                                     rhs = [ pathsK_nodes[i][j][nStations]*float(commodities[i][2]) ] ) # Right Hand Side of the constraint

    try:
        print("\n\n-----------RESTRICTED DUAL SOLUTION : -----------\n")
        model.set_results_stream(None)
        model.set_warning_stream(None)
        model.solve()
        model.write('test_dual.lp')
    except CplexSolverError as e:
        print("Exception raised during dual of restricted problem: " + e)

    # Return the dual variables corresponding to the commodities, arcs and nodes constraints
    return model.solution.get_values()[:nCommodities] , model.solution.get_values()[nCommodities:]

def restricted_Master(pathsK_nodes):
    
#   ------------- SOLVE THE RESTRICTED MASTER PROBLEM -------------
        
    model = cplex.Cplex() # Initialize the model
    model.objective.set_sense(model.objective.sense.minimize) ## Set the objective function to minimization
    
    #Create decision variables
    for i in range(nCommodities):
        for j in range(len(pathsK_nodes[i])):
            model.variables.add(obj = [pathsK_nodes[i][j][nStations]*float(commodities[i][2])], lb = [0], ub = [1],#obj contains the coefficients of the decision variables in the objective function                               
                                types = [model.variables.type.continuous],
                                names = ['P(' + str(i) + ',' + str(j) + ')']) 

    #Add constraints
    #Flow conservation constraints :
    count = 0 # To iterate over the indices of the decision variables
    for i in range(nCommodities):
        ind = [] # Put the indices of the non-zero variables of the current constraint
        val = [] # Put their coefficients here
        row = []
        for j in range(len(pathsK_nodes[i])):
            ind.append(count) 
            val.append(1)
            count += 1
        row.append([ind, val])
        model.linear_constraints.add(lin_expr = row,
                                     senses = "E", #Equality constraint
                                     rhs = [1] ) # Right Hand Side of the constraint
 
    #Node capacity constraints :
    for i in range(nStations): 
        ind, val, row = [], [], []
        count = 0
        for j in range(nCommodities): #For each commodity path, check each time a path contains the node (do not include starting node)
            for k in range(len(pathsK_nodes[j])):
                if(pathsK_nodes[j][k][i] > 0): # If it is the case, add the decision variable index to the constraint
                    ind.append(count)
                    val.append(commodities[j][2]) # With its coefficiant
                count += 1
        row.append([ind, val])
        model.linear_constraints.add(lin_expr = row,
                                     senses = "L", # Less-than
                                     rhs = [ capacity_node[node[i][2]-1] ] ) #Capacity of the node (-1 because type_id start at 1)
        
    try:
        print("\n\n-----------RESTRICTED MASTER SOLUTION : -----------\n")
        model.set_results_stream(None)
        model.set_warning_stream(None)
        model.solve()
        model.write('test_restricted.lp')
        print("\n")
        print("Solution primal : ",model.solution.get_values())
        
        #Print the solution
        count = 0
        for i in range(nCommodities):
            for j in range(len(pathsK_nodes[i])):
                indices = [k for k, x in enumerate(pathsK_nodes[i][j]) if (x > 0 and k < nStations)] #Get the indices of the nodes contained in the current path
                print("\t", model.solution.get_values(count)*commodities[i][2] ,"units of commodity n°", i+1 ,"on path", node[commodities[i][0]][0],''
                      + ' '.join([node[k][0] for k in indices]) + ". Length path : " + str(pathsK_nodes[i][j][nStations]) )
                count += 1
                
        print("\nTotal cost = " + str(model.solution.get_objective_value()))
        
        dualCommodities, dualStations = dual_of_Restrited(pathsK_nodes) # Compute the dual variables and print them
        
        print()
        for i in range(len(dualCommodities)):
            if(dualCommodities[i] != 0):
                print("Dual values y_K" + str(i+1) + " = "+ str(dualCommodities[i]))
        for i in range(len(dualStations)):
            if(dualStations[i] != 0):
                print("Dual values y_Node" + str(i+1) + " = "+ str(dualStations[i]))
            
    except CplexSolverError as e:
        print("Exception raised during restricted master problem: ", e)
        return [-1] * 4 # return -1 to indicate the infeasibility of restricted master problem

    return model.solution.get_objective_value(), model.solution.get_values(), dualCommodities, dualStations

def pricingProblem(dualCommodities, dualStations): 
    """ Return the paths corresponding to the most negative reduced cost of each commodity
    """
#   ------------- SOLVE THE PRICING PROBLEM-------------
    
    reducedCost = 0
    forCommodity = [] # id of the commodity which will receive the new path, if any
    bestPath = [] # Contains the path which will be added (arc_id's and cost)
    for i in range(nCommodities): #Solve the shortest path problem (with updated length) for each commodity
        model = cplex.Cplex()
        model.objective.set_sense(model.objective.sense.minimize) ## Set the objective function to minimization
    
        #Create decision variables (array of nArcs size)
        for j in range(nArcs):
            model.variables.add(obj = [commodities[i][2]*(arc[j][2] - dualStations[ arc[j][1]-1 ] )],
                                lb = [0], ub = [1],
                                types = [model.variables.type.integer],
                                names = ['alpha ( ' + str(j) + ' )'])
            
        model.objective.set_offset(-dualCommodities[i]) # Check if does what we want

        #Add the consistency constraint ( decision varibales form a path )
        #For each node, check if it is a starting node, ending node or in-between node
        for k in range(nStations):
            ind = [] 
            val = [] 
            row = []
            if(commodities[i][0] != k and commodities[i][1] != k ): #If the node is in between 
                rhs = [0]
                for j in range(nArcs):
                    if(arc[j][0]-1 == k):# Compute its leaving flow (-1 because in the data we start with index 1 and not 0)
                        ind.append(j) 
                        val.append(1) 
                    elif(arc[j][1]-1 == k): # Minus what comes in
                        ind.append(j)
                        val.append(-1)
            elif(commodities[i][0] == k): # If the node is the starting node of the commodity
                rhs = [1]
                for j in range(nArcs):
                    if(arc[j][0]-1 == k):# Compute its leaving flow
                        ind.append(j) 
                        val.append(1) 
                    elif(arc[j][1]-1 == k): # Minus what comes in
                        ind.append(j)
                        val.append(-1)
            elif(commodities[i][1] == k): # If the node is the starting node of the commodity
                rhs = [-1]
                for j in range(nArcs):
                    if(arc[j][0]-1 == k):# Compute its leaving flow
                        ind.append(j) 
                        val.append(1) 
                    elif(arc[j][1]-1 == k): # Minus what comes in
                        ind.append(j)
                        val.append(-1)                       
            row.append([ind, val])
            model.linear_constraints.add(lin_expr = row,
                                             senses = "E", #Equality constraint
                                             rhs = rhs )
            
        try:
            print("\n\n-----------PRICING PROBLEM SOLUTION FOR COMMODITY n°", i ,": -----------\n")
            model.set_results_stream(None)
            model.set_warning_stream(None)
            model.solve()
            model.write('test_princing.lp')
            #Print reduced cost for the current commodity
            print()
            print("\n\tREDUCED COST : ", model.solution.get_objective_value())
            #print("\n\tNEW PATH : ", model.solution.get_values())
        except CplexSolverError as e:
            print("Exception raised during pricing problem: " + e)

        if(round(model.solution.get_objective_value()) < 0): # If we obtained a negative reduced cost, take it, again round to avoid computational mistake
            reducedCost = model.solution.get_objective_value()
            toAdd = model.solution.get_values()
            forCommodity.append(i)
    
            tempCost = 0
            for i in range(nArcs):
                if(toAdd[i] == 1):
                    tempCost += arc[i][2]
            toAdd.append(tempCost) #Put the cost of the path at the bottom of its arcs description
            
            bestPath.append(toAdd)
        
    return reducedCost, bestPath, forCommodity
                   
# ------------- DATA ------------- 

s_time_loading_data = time.time() # Time - 0 Start loading data

nStations, node = getNodes() # Get the nodes
nArcs, arc = getArcs() # Get the arcs
nCommodities, commodities = getCommodities() # Get the commodities
cost = getCosts(nStations) # Get the cost adjency matrix

scale_capacity = 1 # Variable used to scale up or down the capacity of the arcs and nodes (depending on their type)
capacity_node = [45, 80, 150, 200] # Define a capacity value for each type of node
capacity_node = [cap*scale_capacity for cap in capacity_node] # Scale them by a constant

t_time_loading_data = time.time() # Time - 1 Start searching for feasible solution

# ------------- ALGORITHM -------------

shortest_paths = np.zeros((nStations, nStations), dtype=int) # Represent the shortest path distance between each nodes
previous_nodes = np.zeros((nStations, nStations), dtype=int) # Represent the previous arc of each node of the shortest path from index

#pathsK_nodes[k,i,j] with k = commodity_id // i = path_id of commodity k // j = node_id of the path 
pathsK_nodes = [[] for i in range(nCommodities)] # Contains all initial paths for each commodity (by nodes)

usedNode = set() # Contains the node deleted from the network
mostUsedNode = -1
pathsK_nodes = getInitSet() # Get the initial set of variable through Dijkstra and "Ignoring most violated node" method

initFind = False # To know if feasible initial set has been found

while True:
    #Check if init set is feasible
    while True:
        obj_Function, solution, dualCommodities, dualStations = restricted_Master(pathsK_nodes) #Solve the restricted master problem
        if(obj_Function > -1):  # We found a feasible solution
            if(initFind == False): # Just to set the timer correctly
                t_time_init_set = time.time() # Time - 2 Start optimizing the problem
                initFind = True
            break
        mostUsedNode = getMostUsedNode() #Compute the most violated node
        usedNode.add(mostUsedNode)
        print('IGNORE : ', mostUsedNode)
        pathsK_nodes= getInitSet()
    
    #Then iterate over pricing and restricted master problem
    reducedCost, newPath, forCommodity = pricingProblem(dualCommodities, dualStations) #Solve the pricing problem
    
    if(round(reducedCost) == 0): # round to avoid minor computational error (10^-23 for 0)
        t_time_solving = time.time() # Time - 3 Algorithm terminates
        #Print final solution
        count = 0
        for i in range(nCommodities):
            print()
            for j in range(len(pathsK_nodes[i])):
                
                indices = [k for k, x in enumerate(pathsK_nodes[i][j]) if (x > 0 and k < nStations) ] #Get the indices of the arcs contained in the current path
                print("\t", solution[count]*commodities[i][2] ,"units of commodity n°", i+1 ,"on path", node[commodities[i][0]][0],''
                      + ' '.join([node[k][0] for k in indices]) + ". Length path : " + str(pathsK_nodes[i][j][nStations]) )
                count += 1
                
        print("\nTotal cost = " + str(obj_Function))
        print("\nLoading data duration : ", t_time_loading_data - s_time_loading_data)
        print("Finding initial set duration : ", t_time_init_set - t_time_loading_data)
        print("Solving problem duration : ", t_time_solving - t_time_init_set)
        print("Total duration : ", t_time_solving - s_time_loading_data)
        break
    else:
        for j in range(len(newPath)):
            newPath_nodes = [0] * nStations #newPath is an array of arcs, transform it into an array of stations
            for i in range(nArcs):
                if(newPath[j][i] == 1):
                    newPath_nodes[arc[i][1]-1] = 1
                    
            newPath_nodes.append(newPath[j][nArcs]) #Append the cost of the path
            pathsK_nodes[forCommodity[j]].append(newPath_nodes) #Finally append the path to initial set
