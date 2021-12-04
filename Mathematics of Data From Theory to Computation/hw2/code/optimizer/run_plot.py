import matplotlib.pyplot as plt
import pickle
import os
os.environ['MKL_SERVICE_FORCE_INTEL'] = '1'
choice_dir = '1e-2/seed_0'
outdir = os.path.join('./output/', choice_dir)
if not os.path.exists(outdir):
    os.makedirs(outdir)

os.system('python main.py --optimizer sgd --learning_rate 1e-2 --output=' +
          outdir + '/sgd.pkl')
os.system(
    'python main.py --optimizer momentumsgd --learning_rate 1e-2 --output=' +
    outdir + '/momentumsgd.pkl')
os.system('python main.py --optimizer rmsprop --learning_rate 1e-2 --output=' +
          outdir + '/rmsprop.pkl')
os.system('python main.py --optimizer amsgrad --learning_rate 1e-2 --output=' +
          outdir + '/amsgrad.pkl')
optimizers = ['sgd', 'momentumsgd', 'rmsprop', 'amsgrad']

fig_dir = os.path.join('./figs/', choice_dir)
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)

# Plots the training losses.
plt.figure()
for optimizer in optimizers:
    data = pickle.load(open(outdir + '/' + optimizer + ".pkl", "rb"))
    plt.plot(data['train_loss'], label=optimizer)
plt.ylabel('Trainig loss')
plt.xlabel('Epochs')
plt.legend()
plt.savefig(os.path.join(fig_dir, 'loss.pdf'))
plt.show()

# Plots the training accuracies.
plt.figure()
for optimizer in optimizers:
    data = pickle.load(open(outdir + '/' + optimizer + ".pkl", "rb"))
    plt.plot(data['train_accuracy'], label=optimizer)
plt.ylabel('Trainig accuracy')
plt.xlabel('Epochs')
plt.legend()
plt.savefig(os.path.join(fig_dir, 'accuracy.pdf'))
plt.show()
