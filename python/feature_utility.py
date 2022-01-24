from math import floor, log10
import numpy as np
from numpy import zeros


def frq2mm(frequency):
    a = 0.06
    k = 165.4
    mm = (1.0 / a) * log10(frequency / k + 1.0)
    return mm


def mm2frq(mm):
    a = 0.06
    first = np.multiply(mm, a)
    frequency = 165.4 * (np.power(10.0, first) - 1.0)
    return frequency


def greenwood(num_filters, lowest, highest):
    base = np.arange(0, num_filters + 1, 1)
    places = base * ((frq2mm(highest) - frq2mm(lowest)) / num_filters)
    places = [base + frq2mm(lowest) for base in places]

    freqs = mm2frq(places)
    lower = freqs[0:num_filters]
    upper = freqs[1:num_filters + 1]

    return lower, upper


def reduce_filter_num(info, num_filters, feature_idx, fs):
    data = info["data"]
    _, num_obs = data.shape
    _, orig_filters = data[feature_idx, 0].shape
    max_freq = floor(fs / 2)
    step_size = max_freq / orig_filters
    lower, upper = greenwood(num_filters, 1, max_freq)

    for obs in range(0, num_obs):
        obs_data = data[feature_idx, obs]
        num_samples, orig_filters = obs_data.shape
        new_mat = zeros((num_samples, num_filters))
        for i in range(0, num_samples):
            sample_data = obs_data[i, :]
            for j in range(0, num_filters):
                bottom_filt = floor(lower[j] / step_size)
                top_filt = floor(upper[j] / step_size)
                if top_filt > orig_filters:
                    top_filt = orig_filters
                vals = sample_data[bottom_filt:top_filt+1]
                total_freq = np.sum(vals)
                new_mat[i, j] = total_freq / ((top_filt - bottom_filt) + 1)
        obs_data = new_mat
        data[feature_idx, obs] = obs_data

    return data

def restore_original_filters(feature, orig_filters, fs):
    num_samples, num_filters = feature.shape
    if num_filters == orig_filters:
        return feature
    max_freq = floor(fs / 2)
    step_size = max_freq / orig_filters
    lower, upper = greenwood(num_filters, 1, max_freq)

    new_mat = zeros((num_samples, orig_filters))
    for i in range(0, num_samples):
        sample_data = feature[i, :]
        for j in range(0, num_filters):
            bottom_filt = floor(lower[j] / step_size)
            top_filt = floor(upper[j] / step_size)
            if top_filt > orig_filters:
                top_filt = orig_filters
            val = sample_data[j]
            for x in range(bottom_filt, top_filt):
                new_mat[i, x] = val
    feature = new_mat

    return feature

