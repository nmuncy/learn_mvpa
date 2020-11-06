from mvpa2.base import cfg
import numpy as np
import pylab as pl

sumo_wrestlers_height = [4, 2, 2, 3, 4]
sumo_wrestlers_weight = [8, 6, 2, 5, 7]
basketball_players_height = [3, 4, 5, 5, 3]
basketball_players_weight = [2, 5, 3, 7, 3]

pl.plot(sumo_wrestlers_height, sumo_wrestlers_weight, 'ro',
        linewidth=2, label="Sumo wrestlers")
pl.plot(basketball_players_height, basketball_players_weight, 'bx',
        linewidth=2, label="Basketball players")
pl.xlim(0, 6)
pl.ylim(0, 10)
pl.xlabel('Height')
pl.ylabel('Weight')
pl.legend()


# transpose to have observations along the first axis
sumo_data = np.vstack((sumo_wrestlers_height,
                       sumo_wrestlers_weight)).T
# same for the baskball data
basketball_data = np.vstack((basketball_players_height,
                             basketball_players_weight)).T
# now stack them all together
all_data = np.vstack((sumo_data, basketball_data))


# creates: [  1,  1,  1,  1,  1, -1, -1, -1, -1, -1]
all_desired_output = np.repeat([1, -1], 5)
zero_meaned_data = all_data - all_data.mean(axis=0)


# compute pseudo-inverse as a matrix
pinv = np.linalg.pinv(np.mat(zero_meaned_data))
# column-vector of observations
y = all_desired_output[np.newaxis].T
weights = pinv * y
gridspec = np.linspace(-4, 4, 20)
input_space_X, input_space_Y = np.meshgrid(gridspec, gridspec)


# for the rest it is easier to have `weights` as a simple array, instead
# of a matrix
weights = weights.A
weighted_output_Z = input_space_X * weights[0] + input_space_Y * weights[1]
pl.figure()
pl.pcolor(input_space_X, input_space_Y, weighted_output_Z,
          cmap=pl.cm.Spectral)
pl.plot(zero_meaned_data[all_desired_output == 1, 0],
        zero_meaned_data[all_desired_output == 1, 1],
        'ro', linewidth=2, label="Sumo wrestlers")
pl.plot(zero_meaned_data[all_desired_output == -1, 0],
        zero_meaned_data[all_desired_output == -1, 1],
        'bx', linewidth=2, label="Basketball players")
pl.xlim(-4, 4)
pl.ylim(-4, 4)
pl.colorbar()
pl.xlabel('Demeaned height')
pl.ylabel('Demeaned weight')
pl.title('Decision output')
pl.legend()
