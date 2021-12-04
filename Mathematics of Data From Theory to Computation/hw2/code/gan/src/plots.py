import matplotlib.pyplot as plt
import imageio

# might need to change the following line depending on your OS
plt.switch_backend('agg')


def compare_samples_2D(sample_1,
                       sample_2,
                       filename,
                       x_lim=[-4, 4],
                       y_lim=[-4, 4]):
    """Plot real vs generated data in 2D"""
    plt.scatter(x=sample_1[:, 0], y=sample_1[:, 1], color='blue', alpha=.1)
    plt.scatter(x=sample_2[:, 0], y=sample_2[:, 1], color='red', alpha=.2)
    axes = plt.gca()
    axes.set_xlim(x_lim)
    axes.set_ylim(y_lim)
    axes.set_aspect('equal')
    if type(filename) is list:
        for fl in filename:
            plt.savefig(fl, bbox_inches='tight')
    else:
        plt.savefig(filename, bbox_inches='tight')
    plt.clf()


def animate(input_files, output_file, duration):
    """Generate an animated gif from a list of files"""
    images = [None] * len(input_files)
    for i, file in enumerate(input_files):
        images[i] = imageio.imread(file)
    imageio.mimsave(output_file, images, 'GIF', duration=duration)
