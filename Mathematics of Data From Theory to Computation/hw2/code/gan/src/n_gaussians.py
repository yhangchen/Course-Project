import torch as pt
import numpy as np
from torch.utils.data import Dataset


class NGaussians(Dataset):
    def __init__(self,
                 N=2,
                 moments=None,
                 inter_distance=4.0,
                 max_cols=None,
                 max_examples=None,
                 scaling=0.5):
        """

        :param N:
        :param moments:
        :param inter_distance:
        :param max_cols:
        :param max_examples: by default dim**2 *100
        """
        if moments is not None:
            assert len(moments) == N
            assert all([len(x) == 2 for x in moments])
            loc = pt.stack([
                x[0] if pt.is_tensor(x[0]) else pt.tensor(x[0])
                for x in moments
            ])
            scale = pt.stack([
                x[1] if pt.is_tensor(x[1]) else pt.tensor(x[1])
                for x in moments
            ])
        else:
            if max_cols is None:
                max_cols = int(np.ceil(np.sqrt(N)))
            x = pt.tensor([inter_distance * (i % max_cols) for i in range(N)])
            y = pt.tensor([inter_distance * (i // max_cols) for i in range(N)])
            loc = pt.stack([x, y], -1)
            scale = pt.ones_like(loc)

        if max_examples is None:
            max_examples = int(10**loc.shape[-1] * loc.shape[0])
        loc = loc * scaling
        scale = scale * scaling

        mix = pt.distributions.Categorical(logits=pt.ones(len(loc)))
        comp = pt.distributions.Independent(
            pt.distributions.Normal(loc, scale), 1)
        self.dist = pt.distributions.MixtureSameFamily(mix, comp)
        self.max_examples = max_examples
        self.examples = self.dist.sample([max_examples])

    def __len__(self):
        return self.max_examples

    def __getitem__(self, item):
        return self.examples[item]

    def sample(self, *args, **kwargs):
        return self.dist.sample(*args, **kwargs)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    ds = NGaussians()
    plt.scatter(ds.examples[:, 0], ds.examples[:, 1])
    plt.show()
    ds = NGaussians(N=16)
    plt.scatter(ds.examples[:, 0], ds.examples[:, 1])
    plt.show()
