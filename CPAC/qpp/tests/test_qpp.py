
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
from CPAC.qpp.qpp import detect_qpp

np.random.seed(10)


def test_qpp():

    voxels, trs = 600, 200
    window_length = 15

    x1 = np.sin(2 * np.pi * 10 * np.linspace(0, 1, trs))
    x = np.tile(x1, (voxels, 1)) + np.random.uniform(0, 1, (voxels, trs))
    x -= x.mean()
    x /= x.std()

    scipy.io.savemat('CPAC/qpp/tests/matlab/qpp.mat', mdict={
        'B': x,
        'msk': np.ones((voxels))
    })

    best_template_segment, best_selected_peaks, best_template_metrics = detect_qpp(
        data=x,
        num_scans=4,
        window_length=window_length,
        permutations=1,
        correlation_threshold=0.3,
        iterations=1,
        convergence_iterations=1
    )

    plt.figure()
    plt.fill_between(range(trs), x.mean(axis=0) - x.std(axis=0), x.mean(axis=0) + x.std(axis=0), facecolor='blue', alpha=0.2)
    plt.plot(range(trs), x.mean(axis=0), label="Data", color='blue')

    convolved = np.zeros((voxels, trs))
    for i in range(trs):
        convolved[i] = (np.convolve(x[i], best_template_segment[i], mode='full') / window_length)[:trs]

    plt.fill_between(range(trs), convolved.mean(axis=0) - convolved.std(axis=0), convolved.mean(axis=0) + convolved.std(axis=0), facecolor='red', alpha=0.2)
    plt.plot(range(trs), convolved.mean(axis=0), label="Convolved QPP", color='red')

    plt.plot(range(-window_length, 0), best_template_segment.mean(axis=0), '--', label="QPP", color='green')
    plt.axvline(x=0, color='black')

    plt.legend()
    plt.show()

    plt.figure()
    plt.fill_between(range(trs), x.mean(axis=0) - x.std(axis=0), x.mean(axis=0) + x.std(axis=0), facecolor='blue', alpha=0.2)
    plt.plot(range(trs), x.mean(axis=0), label="Data", color='blue')
    for xc in best_selected_peaks:
        plt.axvline(x=xc, color='r')
    plt.legend()
    plt.show()