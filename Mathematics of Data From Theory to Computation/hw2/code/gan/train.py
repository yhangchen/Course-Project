import torch
import argparse

from torch import distributions

from src.n_gaussians import NGaussians
from src.trainer import GanTrainer
from src.variables import DualVariable, Generator
from src.optim import Adam

import matplotlib.pyplot as plt


def main(args):
    torch.manual_seed(args.seed)
    device = torch.device(args.device)

    # Define true distribution and noise distrubution
    # feel free to play with these to see the effects
    data = NGaussians(N=2)
    noise_mean = torch.zeros(2)
    noise_covariance = torch.eye(2)
    noise = distributions.MultivariateNormal(noise_mean, noise_covariance)

    # plot the real data
    plt.scatter(data.examples[:, 0], data.examples[:, 1])
    plt.show()

    # Initialize generator and dual variable. Again, feel free to play with the hidden_dim
    f = DualVariable(input_dim=2, hidden_dim=args.hidden_dim, c=args.clip)
    g = Generator(noise_dim=2, output_dim=2, hidden_dim=args.hidden_dim)

    # Create optimizers and initialize adam optimizers with betas=(0.0,0.9), one for generator and dual each
    # For API, see our pytorch tutorial or https://pytorch.org/docs/stable/optim.html
    f_optim = Adam(f.parameters(), lr=args.lr, betas=(0.0, 0.9))
    g_optim = Adam(g.parameters(), lr=args.lr, betas=(0.0, 0.9))

    # Initialize trainer
    trainer = GanTrainer(args.batch_size,
                         data=data,
                         noise=noise,
                         make_gif=args.make_gif)

    # train and save GIF
    trainer.alternating(n_iter=args.n_iter,
                        f=f,
                        g=g,
                        f_optim=f_optim,
                        g_optim=g_optim,
                        f_ratio=args.f_ratio,
                        n_checkpoints=4)
    trainer.render_gif('movie' + '.gif', duration=0.1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # training parameters
    parser.add_argument('--batch_size', type=int, default=200)
    parser.add_argument('--n_iter', type=int, default=3000)
    parser.add_argument('--hidden_dim', type=int, default=100)
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--clip', type=float, default=1e-2)
    parser.add_argument('--f_ratio', type=int, default=5)

    # miscelaneous parameters
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--make_gif', type=bool, default=True)
    parser.add_argument('--device', type=str, default='cpu')
    args = parser.parse_args()
    main(args)
