## Function to add noise
import numpy as np

def Add_noise(data, reporting_rate, phi, seed=None):

    if seed is not None:
        np.random.seed(seed)
    
    data = np.array(data)
    observed_cases = []
    for i in range(data.shape[0]):
        case = []
        for t in range(data.shape[1]):
            mean_cases = reporting_rate * data[i, t]

            # Negative binomial: r = phi (dispersion), p = phi/(phi + mean)
            p = phi / (phi + mean_cases)

            case.append(np.random.negative_binomial(phi, p))
        observed_cases.append(case)
    return np.array(observed_cases)