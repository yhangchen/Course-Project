import torch

from torch import nn


class Generator(nn.Module):
    def __init__(self, noise_dim=2, output_dim=2, hidden_dim=100):
        super().__init__()
        self.inner = nn.Sequential(
            nn.Linear(noise_dim, hidden_dim, bias=True), nn.ReLU(),
            nn.Linear(hidden_dim, output_dim, bias=True))

    def forward(self, z):
        """
        Evaluate on a sample. The variable z contains one sample per row
        """
        return self.inner(z)


class DualVariable(nn.Module):
    def __init__(self, input_dim=2, hidden_dim=100, c=1e-2):
        super().__init__()
        self.c = c
        self.inner = nn.Sequential(nn.Linear(input_dim, hidden_dim, bias=True),
                                   nn.ReLU(),
                                   nn.Linear(hidden_dim, 1, bias=True))

    def forward(self, x):
        """
        Evaluate on a sample. The variable x contains one sample per row
        """
        return self.inner(x)

    def enforce_lipschitz(self):
        """Enforce the 1-Lipschitz condition of the function by doing weight clipping or spectral normalization"""
        self.spectral_normalisation()
        # <= you have to implement this one
        # self.weight_clipping()
        # <= this one is for another year/only for you as a bonus if you want to compare

    def spectral_normalisation(self):
        """
        Perform spectral normalisation, forcing the singular value of the weights to be upper bounded by 1.
        """
        self.inner.apply(SpectralNorm())

    def weight_clipping(self):
        """
        Clip the parameters to $-c,c$. You can access a modules parameters via self.parameters().
        Remember to access the parameters  in-place and outside of the autograd with Tensor.data.
        """
        self.inner.apply(WeightClipper(self.c))


class WeightClipper():
    def __init__(self, c):
        self.c = c

    def __call__(self, module):
        if isinstance(module, nn.Linear):
            w = module.weight.data
            module.weight.data = w.clamp(-self.c, self.c)


class SpectralNorm():
    def __init__(self, n_power_iterations=1, eps=1e-12):
        self.n_power_iterations = n_power_iterations
        self.eps = eps

    def __call__(self, module):
        if isinstance(module, nn.Linear):
            module = spectral_norm(module,
                                   n_power_iterations=self.n_power_iterations,
                                   eps=self.eps)


def spectral_norm(module, n_power_iterations=1, eps=1e-12):
    weight = getattr(module, 'weight', None)
    torch.nn.utils.parametrize.register_parametrization(
        module, 'weight', _SpectralNorm(weight, n_power_iterations, eps))
    return module


class _SpectralNorm(nn.Module):
    def __init__(self, weight, n_power_iterations=1, eps=1e-12):
        super().__init__()
        ndim = weight.ndim
        self.eps = eps
        if ndim > 1:
            self.n_power_iterations = n_power_iterations
            weight_mat = weight.flatten(1)
            h, w = weight_mat.size()
            u = weight_mat.new_empty(h).normal_(0, 1)
            v = weight_mat.new_empty(w).normal_(0, 1)
            self.register_buffer(
                '_u', torch.nn.functional.normalize(u, dim=0, eps=self.eps))
            self.register_buffer(
                '_v', torch.nn.functional.normalize(v, dim=0, eps=self.eps))
            self._power_method(weight_mat, 15)

    @torch.autograd.no_grad()
    def _power_method(self, weight_mat, n_power_iterations):
        assert weight_mat.ndim > 1
        for _ in range(n_power_iterations):
            self._u = torch.nn.functional.normalize(torch.mv(
                weight_mat, self._v),
                                                    dim=0,
                                                    eps=self.eps,
                                                    out=self._u)
            self._v = torch.nn.functional.normalize(torch.mv(
                weight_mat.t(), self._u),
                                                    dim=0,
                                                    eps=self.eps,
                                                    out=self._v)

    def forward(self, weight):
        if weight.ndim == 1:
            return torch.nn.functional.normalize(weight, dim=0, eps=self.eps)
        else:
            weight_mat = weight.flatten(1)
            if self.training:
                self._power_method(weight_mat, self.n_power_iterations)
            u = self._u.clone(memory_format=torch.contiguous_format)
            v = self._v.clone(memory_format=torch.contiguous_format)
            sigma = torch.dot(u, torch.mv(weight_mat, v))
            return weight / sigma