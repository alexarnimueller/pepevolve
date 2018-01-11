"""
# PEPMORPH

Evolutionary algorithm to morph two given sequences into each other with the help of an evolutionary algorithm. This
script is designed to work with peptide sequences. However, any sequence can be provided, given a suitable distance
matrix is referenced.

:Usage:
python pepmorph.py <start> <target> --lambd <int> --sigma <float> --matrixfile <file> --skip <str> --seed <int>

:Parameters:
:param parent: {str} parent string from which you wish to generate children
:param target: {str} target string to which the parent should be morphed to
:param lambd: {int} number of offsprings to generate in between the two sequences
:param sigma: {float} width of the gaussian distribution for tuning the distance to the previous sequence
:param matrixfile: {str} filename of the distance matrix to use
:param skip: {str} letters (AA) to skip when sampling
:param seed: {int} random seed to use when sampling (makes runs reproducible)

:Example:
python pepmorph.py KLLKLLKKLLKLLKLL GLFDIVKKVVGALGSL --lambd 8 --sigma 0.1 --matrixfile grantham.txt --skip CM --seed 42

:Output:
generated sequences with corresponding distances and sigma values written to the file ``restult.txt``
"""

import numpy as np
import argparse


def shift_sigma(sigma):
    """ Randomly sample a float from a normal distribution to shift the given sigma

    :return: {float} shifted sigma value
    """
    a = np.random.random_sample(1)
    b = np.random.random_sample(1)
    shift = sigma * np.sqrt(-2.0 * np.log10(a)) * np.sin(2.0 * np.pi * b)
    return abs(sigma + shift)[0]


def load_matrix(filename):
    """ Load a smiliarity matrix from the given filename. Values need to be tab-separated with AA as column headers.

    :param filename: {str} filename of the matrix file

    :return: matrix with similarity values, array of corresponding amino acids (file header)
    """
    return np.genfromtxt(filename, delimiter='\t', skip_header=True), np.genfromtxt(filename, max_rows=1, dtype=str)


def get_distance(start, target, matrix, aas):
    """ Calcualte the distance between two sequences given a distance matrix

    :param start: {str} sequence 1
    :param target: {str} sequence 2
    :param matrix: {array} distance matrix for all letters
    :param aas: {array} list of letters corresponding to the matrix columns
    :return: {float} distance
    """
    distance = float()
    for i in range(len(start)):
        a = np.where(start[i] == aas)[0]
        b = np.where(target[i] == aas)[0]
        distance += matrix[a, b] ** 2
    return np.sqrt(distance)[0]


def get_step_dist(start, target, lambd, matrix, aas):
    """ Get the number of steps and the distance per step from start to target sequence

    :param start: {str} start string to morph from
    :param target: {str} target string to morph to
    :param lambd: {int} number of steps to take in between
    :param matrix: {array} distance matrix for all letters
    :param aas: {array} list of letters corresponding to the matrix columns
    :return: {array} array of length ``lambd`` with distances from the target sequences
    """
    # get total distance between start and target
    distance = get_distance(start, target, matrix, aas)

    # divide distance by nr of steps
    step_dist = distance / float(lambd + 1.)

    return np.arange(step_dist, distance, step_dist)[::-1].round(1), step_dist


def mutate(parent, sigma, matrix, aas, skip_aa=None):
    """ Mutate a given parent sequence with given sigma and distance matrix to a child through an evolutionary algorithm.

    :param parent: {str} parent sequence.
    :param sigma: {float} sigma used to generate gaussian
    :param matrix: {array} distance matrix to use
    :param aas: {list} array of amino acids corresponding to the columns in ``matrix``
    :param skip_aa: {str} amino acids / characters to skip. If more then one, just append, e.g. ``CM`` for C and M
    :return: {str} generated child and corresponding distance (if ``distance=True``)
    """
    child = str()
    distance = float()
    sigma = shift_sigma(sigma)

    # min-max scale the matrix
    matrix = (matrix - np.min(matrix, axis=1)) / (np.max(matrix, axis=1) - np.min(matrix, axis=1))

    for p in range(len(parent)):
        # get the mutation probability for this letter in the parent sequence depending on the sigma
        indx = np.where(aas == parent[p])[0][0]
        a = (matrix[indx, :] ** 2) / (2 * sigma ** 2)
        b = np.exp(-(matrix[indx, :] ** 2) / (2 * sigma ** 2))
        prob = np.exp(-a / np.sum(b))

        # handle unwanted letters
        aa_flag = np.ones(len(aas))
        if skip_aa:
            for aa in skip_aa:
                idx_aa = int(np.where(aa == aas)[0])
                aa_flag[idx_aa] = 0
        probas = prob * aa_flag / (np.sum(prob * aa_flag))

        # sample letter with given probability
        sample = np.random.choice(aas, 1, p=probas)
        mut_indx = np.where(sample == aas)[0][0]
        child += sample[0]
        distance += matrix[indx, mut_indx] ** 2

    return child, np.sqrt(distance), sigma


def main(start, target, lamb, sig, filename, skip_aa=None):
    matrix, aas = load_matrix(filename)
    dist_list, step_size = get_step_dist(start, target, lamb, matrix, aas)
    children = list()

    with open('result.txt', 'w') as f:
        f.write("PEPMORPH RESULTS\n================\n\n"
                "Sigma:\t%.5f\nLambda:\t%i\nSkip:\t%s\nStart:\t%s\nTarget:\t%s\n\nDist\tSigma\tSequence\n" %
                (sig, lamb, skip_aa, start, target))

        parent = start
        for l in range(lamb):
            print("Running step %i..." % (l + 1))
            child, _, used_sig = mutate(parent, sig, matrix, aas, skip_aa)
            dist = get_distance(target, child, matrix, aas)

            while abs(dist - dist_list[l]) > (step_size / 4.):  # while achieved distance isn't in an acceptable range
                child, _, used_sig = mutate(parent, sig, matrix, aas, skip_aa)
                dist = get_distance(target, child, matrix, aas)

            children.append(child)
            f.write(str(dist.round(3)) + "\t" + str(used_sig.round(3)) + "\t" + child + "\n")
            parent = child


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("start", help="sequence from which to start morphing", type=str)
    parser.add_argument("target", help="target sequence to which to morph", type=str)
    parser.add_argument("-l", "--lambd", help="number of intermediate sequences to produce", type=int, default=10)
    parser.add_argument("-i", "--sigma", help="spread of the gaussian distribution", type=float, default=0.05)
    parser.add_argument("-m", "--matrixfile", help="filename of the distance matrix", type=str, default='grantham.txt')
    parser.add_argument("-n", "--skip", help="letters to skip when sampling", type=str, default='CM')
    parser.add_argument("-s", "--seed", help="random seed to use", type=int, default=42)
    args = parser.parse_args()

    np.random.seed(seed=args.seed)
    if len(args.start) == len(args.target):
        main(start=args.start, target=args.target, lamb=args.lambd, sig=args.sigma,
             filename=args.matrixfile, skip_aa=args.skip)
    else:
        raise IOError("The two provided sequences need to have the same length!")
