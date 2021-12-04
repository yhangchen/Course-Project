import os
import torch
import torch.nn as nn
from torch.autograd import Variable
import torchvision.datasets as dset
import torchvision.transforms as transforms
import torch.nn.functional as F
import torch.optim as optim
from optimizers import createOptimizer
import argparse
import pickle


def main(args):
    torch.manual_seed(args.seed)
    root = './data'
    if not os.path.exists(root):
        os.mkdir(root)

    trans = transforms.Compose([transforms.ToTensor()
                                ])  #, transforms.Normalize((0.5,), (1.0,))])
    # if not exist, download mnist dataset
    train_set = dset.MNIST(root=root,
                           train=True,
                           transform=trans,
                           download=True)
    test_set = dset.MNIST(root=root,
                          train=False,
                          transform=trans,
                          download=True)

    train_loader = torch.utils.data.DataLoader(dataset=train_set,
                                               batch_size=args.batch_size,
                                               shuffle=True)
    test_loader = torch.utils.data.DataLoader(dataset=test_set,
                                              batch_size=args.batch_size,
                                              shuffle=False)

    class MLPNet(nn.Module):
        def __init__(self):
            super(MLPNet, self).__init__()
            self.fc1 = nn.Linear(784, 256)
            self.fc2 = nn.Linear(256, 256)
            self.fc3 = nn.Linear(256, 10)

        def forward(self, x):
            x = x.view(-1, 28 * 28)
            x = F.relu(self.fc1(x))
            x = F.relu(self.fc2(x))
            x = self.fc3(x)
            return x

        def name(self):
            return "MLP"

    ## training
    model = MLPNet()
    optimizer = createOptimizer(args, model)
    criterion = nn.CrossEntropyLoss()
    losses = []
    accuracies = []
    for epoch in range(args.training_epochs):
        # trainning
        avg_cost = 0
        correct_prediction = 0
        total_batch = len(train_loader)
        for batch_idx, (x, target) in enumerate(train_loader):
            optimizer.zero_grad()
            x, target = Variable(x), Variable(target)
            out = model(x)
            loss = criterion(out, target)
            avg_cost += loss.data.numpy() / total_batch
            loss.backward()
            optimizer.step()

            _, pred_label = torch.max(out.data, 1)
            correct_prediction += (pred_label == target.data).sum()
        avg_accuracy = correct_prediction.data.numpy() / (1.0 * total_batch *
                                                          args.batch_size)
        losses.append(avg_cost)
        accuracies.append(avg_accuracy)
        print("Epoch=%d" % (epoch + 1) + " cost={:.4f}".format(avg_cost) +
              " Accuracy={:.4f}".format(avg_accuracy))
    print("Optimization Finished!")

    # Saves the losses and accuracies
    if args.output is not None:
        pickle.dump({"train_loss":losses, "train_accuracy": accuracies},\
                open(args.output, "wb"))

    # Tests the model on the test set.
    correct_prediction = 0
    total_batch = len(test_loader)
    for batch_idx, (x, target) in enumerate(test_loader):
        out = model(x)
        loss = criterion(out, target)
        _, pred_label = torch.max(out.data, 1)
        correct_prediction += (pred_label == target.data).sum()
    print(
        "Test accuracy:",
        correct_prediction.data.numpy() /
        (1.0 * total_batch * args.batch_size))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--optimizer',
                        type=str,
                        default="sgd",
                        help='Name of the optimizer: sgd,\
       momentumsgd, amsgrad, adagrad, rmsprop')
    parser.add_argument('--learning_rate',
                        type=float,
                        default=1e-4,
                        help='Learning rate')
    parser.add_argument('--training_epochs',
                        type=int,
                        default=15,
                        help='Number of training epochs')
    parser.add_argument('--batch_size',
                        type=int,
                        default=100,
                        help='Number of batch sizes')
    parser.add_argument('--delta',
                        type=float,
                        default=1e-8,
                        help='The camping coefficient')
    parser.add_argument('--tau', default=0.9, help='Decaying parameter')
    parser.add_argument('--rho', type=float, default=0.9, help='momentum')
    parser.add_argument('--beta1',
                        type=float,
                        default=0.9,
                        help='first order decaying parameter')
    parser.add_argument('--beta2',
                        type=float,
                        default=0.999,
                        help='second order decaying parameter')
    parser.add_argument('--output',
                        type=str,
                        default=None,
                        help='Output file to save training loss\
       and accuracy.')
    parser.add_argument('--seed', type=int, default=0, help='Seed.')
    args = parser.parse_args()
    main(args)
