import torch
import torchvision
import torch.nn.functional as F

from torchvision import transforms
from torch import nn


class TwoLayerNet(nn.Module):
    def __init__(self, input, hidden, output):
        """
        create two nn.Linear objects and assign them as attributes

        :param input: dimension of the input
        :param hidden: number of hidden neurons
        :param output: dimension of the output
        """
        super().__init__()
        self.linear1 = nn.Linear(input, hidden)
        self.linear2 = nn.Linear(hidden, output)

    def forward(self, x):
        """
        In the forward method we define what is the output of the network
        given an input x. In this example we use the ReLU as our activation function
        """
        x = F.relu(self.linear1(x))
        x = self.linear2(x)
        return x


if __name__ == '__main__':
    net = TwoLayerNet(input=784, hidden=100, output=10)
    x = torch.randn(784)
    result = net(x)
    print('output of the network at input x: ' + str(result))

    train_dataset = torchvision.datasets.MNIST(root='~/data',
                                               train=True,
                                               transform=transforms.ToTensor(),
                                               download=True)

    train_loader = torch.utils.data.DataLoader(dataset=train_dataset,
                                               batch_size=128,
                                               shuffle=True)

    for _, (x, y) in enumerate(train_loader):
        print('batch size: ' + str(x.shape[0]))
        print('input dimension: ' + str(x[0].shape))
        loss_fn = nn.CrossEntropyLoss()
        x = x.view(x.shape[0], -1)  # reshape 28x28 image to a 1x784 vector
        net.zero_grad()  # set the gradients to 0
        output = net(x)
        loss = loss_fn(output, y)
        loss.backward()  # backpropagation
        for p in net.parameters():
            gradient = p.grad
            # perform an update based on the gradient
        break  # stops the for loop. remove this line to iterate through all the data
