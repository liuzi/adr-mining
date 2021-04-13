import torch
import math
import torch.distributions as dist
from pyro.infer import MCMC, NUTS, HMC
# from pyro.infer import MCMC, NUTS, config_enumerate, infer_discrete

torch.cuda.set_device(2)
device = torch.cuda.current_device() 
# print(torch.cuda.current_device())

dtype = torch.float


def model(data):
    coefs_mean = torch.zeros(dim)
    coefs = pyro.sample('beta', dist.Normal(coefs_mean, torch.ones(3)))
    y = pyro.sample('y', dist.Bernoulli(logits=(coefs * data).sum(-1)), obs=labels)
    return y

def run_model():
    true_coefs = torch.tensor([1., 2., 3.])
    data = torch.randn(2000, 3)
    dim = 3
    labels = dist.Bernoulli(logits=(true_coefs * data).sum(-1)).sample()

    hmc_kernel = HMC(model, step_size=0.0855, num_steps=4)
    mcmc_run = MCMC(hmc_kernel, num_samples=500, warmup_steps=100).run(data)
    posterior = EmpiricalMarginal(mcmc_run, 'beta')
    print(posterior.mean)

run_model()

# # Create random input and output data
# x = torch.linspace(-math.pi, math.pi, 2000, device=device, dtype=dtype)
# y = torch.sin(x)

# # Randomly initialize weights
# a = torch.randn((), device=device, dtype=dtype)
# b = torch.randn((), device=device, dtype=dtype)
# c = torch.randn((), device=device, dtype=dtype)
# d = torch.randn((), device=device, dtype=dtype)

# learning_rate = 1e-6
# for t in range(2000):
#     # Forward pass: compute predicted y
#     y_pred = a + b * x + c * x ** 2 + d * x ** 3

#     # Compute and print loss
#     loss = (y_pred - y).pow(2).sum().item()
#     if t % 100 == 99:
#         print(t, loss)

#     # Backprop to compute gradients of a, b, c, d with respect to loss
#     grad_y_pred = 2.0 * (y_pred - y)
#     grad_a = grad_y_pred.sum()
#     grad_b = (grad_y_pred * x).sum()
#     grad_c = (grad_y_pred * x ** 2).sum()
#     grad_d = (grad_y_pred * x ** 3).sum()

#     # Update weights using gradient descent
#     a -= learning_rate * grad_a
#     b -= learning_rate * grad_b
#     c -= learning_rate * grad_c
#     d -= learning_rate * grad_d


# print(f'Result: y = {a.item()} + {b.item()} x + {c.item()} x^2 + {d.item()} x^3')