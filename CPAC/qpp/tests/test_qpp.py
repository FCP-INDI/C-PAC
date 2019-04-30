
import matplotlib.pyplot as plt
import numpy as np

np.random.seed(2)

def test_qpp():
    x1 = np.sin(2 * np.pi * 10 * np.linspace(0, 1, 200))
    x1 -= x1.mean()
    x1 /= x1.std()
    x = np.tile(x1, (10, 1))

    best_template_segment, best_selected_peaks, best_template_metrics = detect_qpp(
        data=x,
        num_scans=1,
        window_length=15,
        permutations=20,
        correlation_threshold=.1,
        max_iterations=100,
        convergence_iterations=1
    )

    plt.figure()
    plt.plot(x[0], label="Data")
    plt.plot(np.convolve(x[0], best_template_segment[0], mode='same'), label="Template")
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(x[0], label="Data")
    for xc in best_selected_peaks:
        plt.axvline(x=xc, color='r')
    plt.legend()
    plt.show()