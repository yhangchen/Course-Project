
                     # ------------- COLUMN GENERATION ALGORITHM FOR THE SMALL INSTANCE ------------- 

from __future__ import print_function

import sys
import time

import cplex
from cplex.exceptions import CplexSolverError
import numpy as np
import openpyxl

from small_data import getNodes, getArcs, getCommodities, getCosts # python script to get the data


def dijkstra(origin):

#   -------- DIJKSTRA'S ALGORITHM --------
    
    labels = np.full((nStations), np.inf) # Contains the cost of the shortest paths from origin with init value = infinity
    labels[origin] = 0
    toTreat = np.array([origin]) # Node that need to be processed
    
    while True:
        currentNode = toTreat[np.argmin(np.take(labels, toTreat))] # Pick the node with lowest shortest path from the set toTreat
        for i in range(nStations):
            if(cost[currentNode, i, 0] > 0): #Among all neighbors of currentNodes
                if(labels[i] > labels[currentNode] + cost[currentNode, i, 0]): # If dual constraint is violated
                    
                    labels[i] = labels[currentNode] + cost[currentNode, i, 0] # Lower it to satisfy the feasibility and complementary slackness
                    toTreat = np.append(toTreat, i) # As we found a better path to go to i, need to check on the neighbors of i

        toTreat = np.delete(toTreat,np.where(toTreat==currentNode)[0]) # Delete the node we just processed
        if(toTreat.size == 0): # If there is no node left to treat, we are sure we obtained the best solution.
            break
    return labels

all_Arc = [] # Will contain all the initials paths found for a given commodity
all_Cost = [] # Will contain the cost of each path in all_Arc

def getAllPath(i,j, optimal): #Get all paths from i to j with length < limit * optimal
    
#   -------- 1/2 GET ALL PATHS FROM i TO j --------
    global all_Arc, all_Cost
    
    all_Arc = []
    all_Cost = []
    visited =[False] * nStations #Keep track of the visited node to avoid cycle
    arc = [0] * nArcs # Write the current path in it (which arcs it contains)
    printAllPath(i, -1 ,j, visited, arc, 0, optimal)
    
    for i in range(len(all_Arc)):
        all_Arc[i].append(all_Cost[i]) # Finally append to each path its corresponding cost

def printAllPath(i, precedent , j, visited, arc, costP, optimal):
    """ Recursive function to find paths from i to j

    Algorithm based on Depth-First-Search : 
    precedent : last node_id we processed
    arc : contains the path in process
    costP : current cost of the path in process
    optimal : shortest path distance between the original i and j

    Basically, explore all possible paths from i, if we reach j without violating the distance criteria, append the path to solution.
    Start at node i, mark node i as visited,
    If i == j : We succesfully reached destination, append the path to all_Arc and its cost to all_Cost

    Else : We recursivey iterate over all the unvisited neighbours k of i, which are not too far, using the adjency matrix (Depth First Search)
           and update the arc(i, k) to 1 if we go to node k

    Finally : Each time we are done with a path, reset arc, costP and visited variables"""
    
#   -------- 2/2 THROUGH A RECURSIVE SEARCH FUNCTION --------
    
    visited[i] = True # Mark the node to avoid cycle
    if(precedent > -1): 
        arc[cost[precedent, i, 1]] = 1 # If you moved from a previous node, add the arc between them to remember the pathway
    
    if(i==j): #If we reach the destination
        all_Cost.append(costP) # Appennd the cost of the succesfull path
        all_Arc.append(arc[:]) # Add the path to the initial set 
        
    else: #Otherwise we iterate over all neighbors unvisited nodes (by depth), only if the future distance is not too far from the shortest path 
        for k in range(nStations):
            if(cost[i, k, 0] > 0 and visited[k]==False and (costP+cost[i, k, 0]) <= limit*optimal) :
                printAllPath(k, i ,j ,visited, arc ,costP + cost[i, k, 0], optimal)
            
    if(precedent > -1): # Once you reached a leaf (no where to go), reset arc and visited variable to their previous value
        arc[cost[precedent, i, 1]] = 0
    visited[i] = False

def getInitSolution_Dijkstra():

#   ------------- GET INITIAL SET OF PROMISSING VARIABLES -------------

    #pathsK[k,i,j] with k = commodity_id // i = path_id of commodity k // j = arc_id of the path 
    pathsK = [ [[]] for i in range(nCommodities)] # Contains all initial paths for each commodity

    
    #Now for each commodity, obtain all possible paths with length < limit * optimal
    for i in range(nCommodities):
        getAllPath(commodities[i][0], commodities[i][1], shortest_paths[commodities[i][0], commodities[i][1]] )
        pathsK[i] = all_Arc # Put all the output in pathsK variable
        
    return pathsK

def dual_of_Restrited(pathsK):
    """ Goal : Obtain the dual variables of the restricted master problem
    """
#   ------------- DUAL OF THE RESTRICTED MASTER PROBLEM -------------
    
    model = cplex.Cplex() # Initialize the model
    model.objective.set_sense(model.objective.sense.maximize) ## Set the objective function to a maximization

    # Add variables
    # Correspond to the flow conservation constraints
    for i in range(nCommodities):
        model.variables.add(obj = [1], lb = [-cplex.infinity], ub = [cplex.infinity],#obj contains the coefficients of the decision variables in the objective function
                            types = [model.variables.type.continuous],
                            names = ['Y( K_' + str(i) + ')'])
        
    # Correspond to the arc capacity constraints
    for i in range(nArcs):
        model.variables.add(obj = [ capacity_arc[arc[i][3]-1] ], lb = [-cplex.infinity], ub = [0],
                            types = [model.variables.type.continuous],
                            names = ['Y( A_' + str(arc[i][0]-1) + '_' + str(arc[i][1]-1) + ')'])

    # Correspond to the node capacity constraints       
    for i in range(nStations):
        model.variables.add(obj = [ capacity_node[node[i][2]-1] ], lb = [-cplex.infinity], ub = [0],
                            types = [model.variables.type.continuous],
                            names = ['Y( ' + str(node[i][0])  + ')'])
        

    #Add constraints
    for i in range(nCommodities):
        for j in range(len(pathsK[i])):
            ind = []# ind : Contains the indices of the non-zero variables of the current constraint
            val = []# val : Contains their coefficients 
            row = []
            ind.append(i)
            val.append(1)
            for k in range(nArcs): # Check which arcs are included in the path
                if(pathsK[i][j][k] == 1):
                    ind.append(nCommodities+k) # Indices of arc variables start after the commodities ones
                    val.append(commodities[i][2])
            for k in range(nStations):# Check which stations are included in the path (do not include starting node)
                for l in range(nArcs):
                    if(pathsK[i][j][l] == 1 and arc[l][1]-1 == k):
                        ind.append(nCommodities+nArcs+k)
                        val.append(commodities[i][2])
            row.append([ind, val])
            model.linear_constraints.add(lin_expr = row,
                                     senses = "L", #Less or equal than constraint
                                     rhs = [ pathsK[i][j][nArcs]*float(commodities[i][2]) ] ) # Right Hand Side of the constraint
            
    try:
        print("\n\n-----------RESTRICTED DUAL SOLUTION : -----------\n")
        model.set_results_stream(None)
        model.set_warning_stream(None)
        model.solve()
        model.write('test_dual.lp')
    except CplexSolverError as e:
        print("Exception raised during dual of restricted problem: " + e)

    # Return the dual variables corresponding to the commodities, arcs and nodes constraints
    return model.solution.get_values()[:nCommodities], model.solution.get_values()[nCommodities:nCommodities + nArcs], model.solution.get_values()[nArcs + nCommodities:]

def restricted_Master(pathsK):
    
#   ------------- SOLVE THE RESTRICTED MASTER PROBLEM -------------
        
    model = cplex.Cplex() # Initialize the model
    model.objective.set_sense(model.objective.sense.minimize) ## Set the objective function to a minimization
    
    #Create decision variables
    for i in range(nCommodities):
        for j in range(len(pathsK[i])):
            model.variables.add(obj = [pathsK[i][j][nArcs]*float(commodities[i][2])], lb = [0], ub = [1],#obj contains the coefficients of the decision variables in the objective function                               
                                types = [model.variables.type.continuous],
                                names = ['P(' + str(i) + ',' + str(j) + ')']) 

    #Add constraints
    #Flow conservation constraints :
    count = 0 # To iterate over the indices of the decision variables
    for i in range(nCommodities):
        ind = [] # Put the indices of the non-zero variables of the current constraint
        val = [] # Put their coefficients here
        row = []
        for j in range(len(pathsK[i])):
            ind.append(count) 
            val.append(1)
            count += 1
        row.append([ind, val])
        model.linear_constraints.add(lin_expr = row,
                                     senses = "E", #Equality constraint
                                     rhs = [1] ) # Right Hand Side of the constraint
        
    #Arc capacity constraints :
    for i in range(nArcs): 
        ind, val, row = [], [], []
        count = 0
        for j in range(nCommodities): #For each commodity path, check each time a path contains the arc, 
            for k in range(len(pathsK[j])):
                if(pathsK[j][k][i] == 1): # If it is the case, add the decision variable index to the constraint
                    ind.append(count)
                    val.append(commodities[j][2]) # With its coefficiant
                count += 1
        row.append([ind, val])
        model.linear_constraints.add(lin_expr = row,
                                     senses = "L", # Less-than
                                     rhs = [ capacity_arc[arc[i][3]-1] ] ) #Capacity of the arc (-1 because type_id start at 1)

    #Node capacity constraints :
    for i in range(nStations): 
        ind, val, row = [], [], []
        count = 0
        for j in range(nCommodities): #For each commodity path, which nodes it contains (do not include starting node)
            for k in range(len(pathsK[j])):
                for l in range(nArcs):
                    if(pathsK[j][k][l] == 1 and arc[l][1]-1 == i): # If it is the case, add the decision variable index to the constraint
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
            for j in range(len(pathsK[i])):
                indices = [k for k, x in enumerate(pathsK[i][j]) if x == 1] #Get the indices of the arcs contained in the current path
                print("\t", model.solution.get_values(count)*commodities[i][2] ,"quantity of commodity n°", i ,"on path", node[commodities[i][0]][0],''
                      + ' '.join([node[arc[k][1]-1][0] for k in indices]) + ". Length path : " + str(pathsK[i][j][nArcs]) )
                count += 1
                
        print("\nTotal cost = " + str(model.solution.get_objective_value()))
        
        dualCommodities, dualArcs, dualStations = dual_of_Restrited(pathsK) # Compute the dual variables and print the non-zero ones
        
        print()
        for i in range(len(dualCommodities)):
            if(dualCommodities[i] != 0):
                print("Dual values y_K" + str(i+1) + " = "+ str(dualCommodities[i]))
        for i in range(len(dualArcs)):
            if(dualArcs[i] != 0):
                print("Dual values y_Arc" + str(i+1) + " = "+ str(dualArcs[i]))
        for i in range(len(dualStations)):
            if(dualStations[i] != 0):
                print("Dual values y_Node" + str(i+1) + " = "+ str(dualStations[i]))
            
    except CplexSolverError as e:
        print("Exception raised during restricted master problem: ", e)
        return [-1] * 5 # return -1 to indicate the infeasibility of restricted master problem

    return model.solution.get_objective_value(), model.solution.get_values(), dualCommodities, dualArcs, dualStations

def pricingProblem(dualCommodities, dualArcs, dualStations):
    """ Return the path corresponding to the most violated constraint
    """
#   ------------- SOLVE THE PRICING PROBLEM-------------
    
    reducedCost = 0
    forCommodity = 0 # id of the commodity which will receive the new path, if any
    bestPath = [] # Contains the path which will be added (arc_id's and cost)
    for i in range(nCommodities): #For each commodity, solve the shortest path problem (with updated cost vector) 
        model = cplex.Cplex()
        model.objective.set_sense(model.objective.sense.minimize) ## Set the objective function to minimization
    
        #Create decision variables (array of nArcs size)
        for j in range(nArcs):
            model.variables.add(obj = [commodities[i][2]*(arc[j][2] - dualArcs[j] - dualStations[ arc[j][1]-1 ] )],
                                lb = [0], ub = [1],
                                types = [model.variables.type.integer],
                                names = ['alpha ( ' + str(j) + ' )'])
        
        model.objective.set_offset(-dualCommodities[i]) 

        #Add the consistency constraint ( decision variables form a path )
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
            #Print reduced cost and potential new path for the current commodity
            print()
            print("\n\tREDUCED COST : ", model.solution.get_objective_value())
            print("\n\tNEW PATH : ", model.solution.get_values())
        except CplexSolverError as e:
            print("Exception raised during pricing problem: " + e)
        else:
            if(round(model.solution.get_objective_value()) < reducedCost): # If we obtained a more violated constraint, take it
                reducedCost = model.solution.get_objective_value()
                bestPath = model.solution.get_values()
                forCommodity = i
        
    if(reducedCost < 0): #If reduced cost negative, compute the cost of the new path
        tempCost = 0
        for i in range(nArcs):
            if(bestPath[i] == 1):
                tempCost += arc[i][2]
        bestPath.append(tempCost) #Put the cost of the path at the bottom of its arcs description
       
    return reducedCost, bestPath, forCommodity
        
def getInitSet():
    """ Launch InitSolution process
    """
    return getInitSolution_Dijkstra()           

# ------------- DATA ------------- 

s_time_loading_data = time.time() # Time - 0 Start loading data

nStations, node = getNodes() # Get the nodes
nArcs, arc = getArcs() # Get the arcs
nCommodities, commodities = getCommodities() # Get the commodities
cost = getCosts(nStations) # Get the cost adjency matrix

scale_capacity = 2.2 # Variable used to scale up or down the capacity of the arcs and nodes (depending on their type)
capacity_arc = [200,350,500,650] # Define a capacity value for each type of arc
capacity_arc = [cap*scale_capacity for cap in capacity_arc] # Scale them by a constant

capacity_node = [1000,2000,3000] # Define a capacity value for each type of node
capacity_node = [cap*scale_capacity for cap in capacity_node] # Scale them by a constant

limit = 1 # Include all paths 'p' satisfying dist(p) < limit*shortest_path

t_time_loading_data = time.time() # Time - 1 Start searching for feasible solution

# ------------- ALGORITHM ------------- 

shortest_paths = np.zeros((nStations, nStations)) # Represent the shortest path distance between each nodes

#Run Dijkstra's Algo for each node to get the shortest path distance
for i in range(nStations-1): #Don't need the last node, if you have already find all the shortest path from every node except the last, you already computed the info for the last node
    shortest_paths[i] = dijkstra(i)
    shortest_paths[:, i] = shortest_paths[i] # Fill the matrix row-wise and column-wise, as the cost is symetric
        
pathsK = getInitSet() # Get the initial set of variable through Dijkstra and "< limit*shortest_path " method
initFound = False #Variable used to know if a feasible init set has been found

while True:
    #Check if init set is feasible
    while True:
        obj_Function, solution, dualCommodities, dualArcs, dualStations = restricted_Master(pathsK) #Solve the restricted master problem
        if(obj_Function > -1):  # We found a feasible solution
            if(initFound == False):
                initFound = True
                t_time_init_set = time.time() # Time - 2 Start optimizing the problem
            break
        limit *= 1.05 # Increase by 5%
        if(limit > 4): # Assume there is no feasible solution if limit > 4
            t_time_solving = time.time()
            print("Total duration : ", t_time_solving - s_time_loading_data)
            sys.exit("NO FEASIBLE SOLUTION")
        pathsK = getInitSet()
    
    #Then iterate over pricing and restricted master problem
    reducedCost, newPath, forCommodity = pricingProblem(dualCommodities, dualArcs, dualStations) #Solve the pricing problem
    
    if(round(reducedCost) == 0): # round to avoid minor computational error (10^-23 instead of 0)
        t_time_solving = time.time() # Time - 3 Algorithm terminates
        #Print final solution
        count = 0
        for i in range(nCommodities):
            print()
            for j in range(len(pathsK[i])):
                
                indices = [k for k, x in enumerate(pathsK[i][j]) if x == 1] #Get the indices of the arcs contained in the current path
                print("\t", solution[count]*commodities[i][2] ,"quantity of commodity n°", i+1 ,"on path", node[commodities[i][0]][0],''
                      + ' '.join([node[arc[k][1]-1][0] for k in indices]) + ". Length path : " + str(pathsK[i][j][nArcs]) )
                count += 1
                
        print("\nTotal cost = " + str(obj_Function))
        print("\nLoading data duration : ", t_time_loading_data - s_time_loading_data)
        print("Finding initial set duration : ", t_time_init_set - t_time_loading_data)
        print("Solving problem duration : ", t_time_solving - t_time_init_set)
        print("Total duration : ", t_time_solving - s_time_loading_data)
        print("\nScaling coefficiant : ", limit)
        print("Nodes capacity : ", capacity_node)
        print("Arcs capacity : ", capacity_arc)
        
        break
    else:
        pathsK[forCommodity].append(newPath) # We iterate with a new path for commodity n° "forCommodity"



