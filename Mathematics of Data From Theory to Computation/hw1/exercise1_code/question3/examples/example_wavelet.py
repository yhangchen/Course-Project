import matplotlib.pyplot as plt
import os, sys
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(CURRENT_DIR))
from common.utils import load_image
from common.operators import Representation_Operator


def main():
    # Load the image
    im_shape = (512, 512)
    im = load_image('data/lauterbrunnen.jpg', im_shape)

    # Wavelet Transform operator
    r = Representation_Operator(m=im_shape[0])

    i_wav = r.W(im).reshape(im_shape)
    i_recon = r.WT(i_wav).reshape(im_shape)

    fig, axs = plt.subplots(1, 3, figsize=(12, 4))
    axs[0].imshow(im, cmap='gray')
    axs[0].set_title('Original')
    axs[1].imshow(abs(i_wav)**0.05, cmap='gray')
    axs[1].set_title('Wavelet coefficients')
    axs[2].imshow(i_recon, cmap='gray')
    axs[2].set_title('Inverse Wavelet transform')

    for ax in axs.flatten():
        ax.set_axis_off()
    plt.tight_layout()
    plt.savefig("examples/wavelet.png")
    plt.show()


if __name__ == "__main__":
    main()
