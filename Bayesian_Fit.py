import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner

# Define the linear model
def linear_model(parameters, x_values):
    slope, intercept = parameters
    return slope * x_values + intercept

# Define the log-likelihood function
def log_likelihood(parameters, x_values, y_values, x_errors, y_errors):
    slope, intercept = parameters
    model_values = linear_model(parameters, x_values)
    total_error = np.sqrt(y_errors**2 + (slope * x_errors)**2)
    sigma_squared = total_error**2
    return -0.5 * np.sum((y_values - model_values)**2 / sigma_squared + np.log(sigma_squared))

# Define the log-prior function
def log_prior(parameters):
    slope, intercept = parameters
    if 0.0 < slope < 10.0 and -200.0 < intercept < 200.0:
        return 0.0
    return -np.inf

# Define the log-posterior function
def log_posterior(parameters, x_values, y_values, x_errors, y_errors):
    log_prior_val = log_prior(parameters)
    if not np.isfinite(log_prior_val):
        return -np.inf
    return log_prior_val + log_likelihood(parameters, x_values, y_values, x_errors, y_errors)

# Input Data
x_values = np.array([1858.5, 3254.19, 240.56, 389.563, 509.978, 1034.54, 1258.27, 1795.13, 519.299, 829.71, 1446.98, 2225.88, 2758.24, 3052.86])
x_errors = np.array([1.1018, 0.4202, 0.1258, 0.0161, 0.0178, 0.2888, 0.3806, 0.6457, 0.0773, 0.1399, 0.5624, 11.8172, 1.4090, 1.9001]) * 1
y_values = np.array([1460.820, 2614.511, 121.7817, 244.6974, 344.2785, 778.9045, 964.057, 1408.013, 351.9320, 609.321, 1120.294, 1764.491, 2204.10, 2447.69])
y_errors = np.array([0.005, 0.010, 0.003, 0.008, 0.012, 0.024, 0.05, 0.03, 0.0021, 0.007, 0.006, 0.014, 0.04, 0.03]) * 1
# x_values = np.array([278.9970, 547.6550, 2813.9900, 3597.0200, 3735.8900, 6128.4900])
# y_values = np.array([1.1572, 1.3567, 1.3881, 1.8728, 1.6225, 1.8517])
# x_errors = np.array([0, 0, 0, 0, 0, 0])
# y_errors = np.array([0.0396, 0.0715, 0.0347, 0.0948, 0.0700, 0.0778])

# Set up the MCMC sampler
num_walkers = 100
num_dimensions = 2
sampler = emcee.EnsembleSampler(num_walkers, num_dimensions, log_posterior, args=(x_values, y_values, x_errors, y_errors))

# Initialize the walkers
initial_positions = np.zeros((num_walkers, num_dimensions))
initial_positions[:, 0] =  np.random.rand(num_walkers) * 10    # initial slope between 1 and 10
initial_positions[:, 1] = -150 + np.random.rand(num_walkers) * 350  # initial intercept between -150 and 200
# these ranges are just for the initial positions of the walkers in the MCMC sampler. The walkers are free to explore beyond these ranges during the sampling process, constrained only by the prior distribution defined in the log_prior function.

# Run the MCMC sampling
num_steps = 8000
sampler.run_mcmc(initial_positions, num_steps, progress=True)

# Get the chain and discard burn-in
chain = sampler.get_chain(discard=500, flat=True) # combining the chains from all walkers into a single chain.

# Plot the colorful 2D corner plot for parameter posterior distributions
labels = ["Slope", "Intercept"]
corner_plot = corner.corner(chain, labels=labels, show_titles=True, title_fmt=".8f", title_kwargs={"fontsize": 12})
corner_plot.savefig("Linear_corner_plot.png")

# Plot the central mean/prediction with transparent uncertainty band
fig, ax = plt.subplots()
ax.errorbar(x_values, y_values, xerr=x_errors, yerr=y_errors, color="black", fmt="s", label="Data")

x_range = np.linspace(0, 4000, 400) # np.linspace(start, stop, num)
y_values_for_each_params = [linear_model(params, x_range) for params in chain] # predictions of our calibrated model at 1000 x points for each set of parameters in the chain
y_values_for_each_params = np.array(y_values_for_each_params)  # Convert list to numpy array
print(y_values_for_each_params.shape) # (total number of samples, number of x points)


# Calculate the median, upper and lower percentiles
median = np.percentile(y_values_for_each_params, 50, axis = 0)
upper = np.percentile(y_values_for_each_params, 97.5, axis = 0)
lower = np.percentile(y_values_for_each_params, 2.5, axis = 0)

# Plot the median, upper and lower percentiles with band
ax.fill_between(x_range, lower, upper, color="blue", alpha=0.2, label='95% Credible Interval')
ax.plot(x_range, median, color="blue", alpha=1, linewidth=2, label='Prediction Mean')
ax.set_xlabel("Channel")
ax.set_ylabel("Energy (keV)")
ax.legend(loc="upper left", fontsize=14)
plt.savefig("Linear_prediction.png")
plt.show()

