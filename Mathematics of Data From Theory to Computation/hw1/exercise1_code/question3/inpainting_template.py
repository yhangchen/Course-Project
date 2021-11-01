import time
import numpy as np
import matplotlib.pyplot as plt
from skimage.restoration import denoise_tv_chambolle
from skimage.metrics import structural_similarity as ssim

from common.utils import *
from common.operators import TV_norm, Representation_Operator, p_omega, p_omega_t, l1_prox


def ISTA(fx, gx, gradf, proxg, params):
    ## TO BE FILLED ##
    """
    Function:  [x, info] = ista(fx, gx, gradf, proxg, parameter)
    Purpose:   Implementation of ISTA.
    Parameter: x0         - Initial estimate.
               maxit      - Maximum number of iterations.
               prox_Lips  - Lipschitz constant for gradient.
               lambda     - regularization factor in F(x)=f(x)+lambda*g(x).
    :param fx:
    :param gx:
    :param gradf:
    :param proxg:
    :param parameter:
    :return:
    """

    method_name = 'ISTA'
    print_start_message(method_name)
    tic_start = time.time()

    # Initialize parameters.
    #### YOUR CODE GOES HERE
    maxit = params['maxit']
    x0 = params['x0']
    lmbd = params['lambda']
    alpha = 1 / params['prox_Lips']

    info = {'itertime': np.zeros(maxit), 'fx': np.zeros(maxit), 'iter': maxit}

    x_k = x0
    for k in range(maxit):
        tic = time.time()

        # Update the iterate
        #### YOUR CODE GOES HERE
        x_next = proxg(x_k - alpha * gradf(x_k), alpha * lmbd)
        if np.linalg.norm(x_next - x_k) / max(
                np.linalg.norm(), 1e-8) < params['stopping_criterion']:
            break
        x_k = x_next

        # Compute error and save data to be plotted later on.
        info['itertime'][k] = time.time() - tic
        info['fx'][k] = fx(x_k) + lmbd * gx(x_k)
        if k % params['iter_print'] == 0:
            print_progress(k, maxit, info['fx'][k], fx(x_k), gx(x_k))

    print_end_message(method_name, time.time() - tic_start)
    return x_k, info


def FISTA(fx, gx, gradf, proxg, params, verbose=False):
    ## TO BE FILLED ##
    """
    Function:  [x, info] = fista(fx, gx, gradf, proxg, parameter)
    Purpose:   Implementation of FISTA (with optional restart).
    Parameter: x0            - Initial estimate.
               maxit         - Maximum number of iterations.
               prox_Lips     - Lipschitz constant for gradient.
               lambda        - regularization factor in F(x)=f(x)+lambda*g(x).
               restart_fista - enable restart.
    :param fx:
    :param gx:
    :param gradf:
    :param proxg:
    :param parameter:
    :return:
    """
    def gradient_scheme_restart_condition(x_k, x_k_next, y_k):
        """
        Whether to restart
        """
        return np.sum((y_k - x_k_next) * (x_k_next - x_k)) > 0

    if verbose:
        method_name = 'FISTAR'
    else:
        method_name = 'FISTA'
    print_start_message(method_name)
    tic_start = time.time()

    # Initialize parameters
    #### YOUR CODE GOES HERE
    maxit = params['maxit']
    x0 = params['x0']
    lmbd = params['lambda']
    alpha = 1 / params['prox_Lips']
    t = 1
    y_k = x0.copy()

    info = {'itertime': np.zeros(maxit), 'fx': np.zeros(maxit), 'iter': maxit}

    x_k = x0
    for k in range(maxit):
        tic = time.time()

        # Update iterate
        #### YOUR CODE GOES HERE
        x_next = proxg(y_k - alpha * gradf(y_k), alpha * lmbd)
        if np.linalg.norm(x_next - x_k) / max(
                np.linalg.norm(x_k), 1e-6) < params['stopping_criterion']:
            break
        if verbose and gradient_scheme_restart_condition(x_k, x_next, y_k):
            y_k = x_k
            t = 1
            continue
        t_next = (1 + np.sqrt(1 + 4 * t**2)) / 2
        y_k = x_next + (t - 1) / t_next * (x_next - x_k)
        x_k = x_next
        t = t_next

        # Compute error and save data to be plotted later on.
        info['itertime'][k] = time.time() - tic
        info['fx'][k] = fx(x_k) + lmbd * gx(x_k)
        if k % params['iter_print'] == 0:
            print_progress(k, maxit, info['fx'][k], fx(x_k), gx(x_k))

    print_end_message(method_name, time.time() - tic_start)
    return x_k, info


def reconstructL1(image, indices, optimizer, params):
    # Wavelet operator
    r = Representation_Operator(m=params["m"])
    indices = params['indices']
    image_vec = image.reshape((params['N'], 1))

    # Define the overall operator
    forward_operator = lambda x: p_omega(r.WT(x), indices)
    adjoint_operator = lambda x: r.W(p_omega_t(x, indices, params['m']))

    # Generate measurements
    b = image_vec[indices]  ## TO BE FILLED ##

    fx = lambda x: np.linalg.norm(b - forward_operator(x))**2 / 2
    gx = lambda x: params['lambda'] * np.linalg.norm(x, 1)
    proxg = lambda x, y: l1_prox(x, params['lambda'] * y)
    gradf = lambda x: adjoint_operator(forward_operator(x) - b)
    x, info = optimizer(fx, gx, gradf, proxg, params)
    return r.WT(x).reshape((params['m'], params['m'])), info


def reconstructTV(image, indices, optimizer, params):
    """
        image: undersampled image (mxm) to be reconstructed
        indices: indices of the undersampled locations
        optimizer: method of reconstruction (FISTA/ISTA function handle)
        params:
    """
    indices = params['indices']
    image_vec = image.reshape((params['N'], 1))
    # Define the overall operator
    forward_operator = lambda x: p_omega(x, indices)
    adjoint_operator = lambda x: p_omega_t(x, indices, params['m']).reshape(
        (params['N'], 1))

    # Generate measurements
    b = image_vec[indices]  ## TO BE FILLED ##

    fx = lambda x: np.linalg.norm(b - forward_operator(x))**2 / 2
    gx = lambda x: params['lambda'] * TV_norm(x)  ## TO BE FILLED ##
    proxg = lambda x, y: denoise_tv_chambolle(x.reshape(
        (params['m'], params['m'])),
                                              weight=params["lambda"] * y,
                                              eps=1e-5,
                                              n_iter_max=50).reshape(
                                                  (params['N'], 1))
    gradf = lambda x: adjoint_operator(forward_operator(x) - b)

    x, info = optimizer(fx, gx, gradf, proxg, params)
    return x.reshape((params['m'], params['m'])), info


# %%

if __name__ == "__main__":

    ##############################
    # Load image and sample mask #
    ##############################
    shape = (256, 256)
    params = {
        'maxit': 200,
        'tol': 10e-15,
        'prox_Lips': 1,  ## TO BE FILLED ##,
        'lambda': None,  ## TO BE FILLED ##,
        'x0': np.zeros((shape[0] * shape[1], 1)),
        'restart_criterion': 'default',  ## TO BE FILLED ##, gradient_scheme,
        'stopping_criterion': 1e-8,  ## TO BE FILLED ##,
        'iter_print': 50,
        'shape': shape,
        'restart_param': 50,
        'verbose': True,
        'm': shape[0],
        'rate': 0.4,
        'N': shape[0] * shape[1]
    }
    PATH = './data/oeschinensee.jpg'
    image = load_image(PATH, params['shape'])

    im_us, mask = apply_random_mask(image, params['rate'])
    indices = np.nonzero(mask.flatten())[0]
    params['indices'] = indices

    # Choose optimization parameters
    #######################################
    # Reconstruction with L1 and TV norms #
    #######################################
    def hyper_sweep(lambda_):
        params['lambda'] = lambda_
        t_start = time.time()
        reconstruction_l1 = reconstructL1(image, indices, FISTA, params)[0]
        t_l1 = time.time() - t_start
        psnr_l1 = psnr(image, reconstruction_l1)
        ssim_l1 = ssim(image, reconstruction_l1)

        t_start = time.time()
        reconstruction_tv = reconstructTV(image, indices, FISTA, params)[0]
        t_tv = time.time() - t_start

        psnr_tv = psnr(image, reconstruction_tv)
        ssim_tv = ssim(image, reconstruction_tv)

        # Plot the reconstructed image alongside the original image and PSNR
        fig, ax = plt.subplots(1, 4, figsize=(20, 6))
        ax[0].imshow(image, cmap='gray')
        ax[0].set_title('Original')
        ax[1].imshow(im_us, cmap='gray')
        ax[1].set_title('Original with missing pixels')
        ax[2].imshow(reconstruction_l1, cmap="gray")
        ax[2].set_title(
            'L1 - PSNR = {:.2f}\n SSIM  = {:.2f} - Time: {:.2f}s'.format(
                psnr_l1, ssim_l1, t_l1))
        ax[3].imshow(reconstruction_tv, cmap="gray")
        ax[3].set_title(
            'TV - PSNR = {:.2f}\n SSIM  = {:.2f}  - Time: {:.2f}s'.format(
                psnr_tv, ssim_tv, t_tv))
        [axi.set_axis_off() for axi in ax.flatten()]
        plt.tight_layout()
        plt.savefig(f'results/log_lambda_{lambda_:3f}.png')  # f-string.
        plt.show()

        out = {}
        out['PSNR'] = (psnr_l1, psnr_tv)
        return out

    log_lambda_list = np.linspace(-4, 1, 21, endpoint=True)
    lambda_list = 10**log_lambda_list
    psnr_l1_list = []
    psnr_tv_list = []
    for log_lambda_ in log_lambda_list:
        lambda_ = 10**log_lambda_
        (psnr_l1, psnr_tv) = hyper_sweep(lambda_)['PSNR']
        psnr_l1_list.append(psnr_l1)
        psnr_tv_list.append(psnr_tv)
    # print(psnr_l1_list)
    # print(psnr_tv_list)
    # psnr_l1_list = [
    #     7.429952778626254, 7.430848138382753, 7.433679107264628,
    #     7.442625609099876, 7.4708658302368285, 7.559757555096498,
    #     7.837229450876196, 8.690528424247201, 11.337757215895856,
    #     19.004290226536856, 21.022878953311537, 20.847370664997694,
    #     20.9381263781152, 21.040782520804264, 20.2792960612071,
    #     17.771801665309344, 13.750049079954376, 6.946897921685781,
    #     5.236148833790418, 5.236148833790418, 5.236148833790418
    # ]
    # psnr_tv_list = [
    #     7.4305396602479625, 7.432704273050945, 7.439551190879086,
    #     7.4612211227924465, 7.529926896043955, 7.748973745227056,
    #     8.458601929049147, 10.84680995246671, 19.00476753023856,
    #     23.159751133430216, 23.546640239895183, 23.541957842801963,
    #     23.3933090810164, 22.623772921815238, 20.625059229600893,
    #     18.575092297243895, 18.043954711382636, 17.871629108843784,
    #     15.89709686106492, 9.48745172035187, -0.0016361003998855198
    # ]
    fig, ax = plt.subplots()
    ax.semilogx(lambda_list, psnr_l1_list, label="L1 reg")
    ax.semilogx(lambda_list, psnr_tv_list, label="TV reg")
    ax.set_xlabel('lambda')
    ax.set_ylabel('PSNR')
    plt.title('PSNR-lambda')
    plt.savefig('results/PSNR-lambda.png')
    plt.show()
