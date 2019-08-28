import numpy as np
import rdp
import pytest
import polyprox


@pytest.fixture()
def random_walk():
    ## Random walk
    N = 2000

    phi = np.random.uniform(0.0, np.pi, N)
    s = np.random.uniform(1.0, 5.0, N)

    x = np.cumsum(s * np.cos(phi))
    y = np.cumsum(s * np.sin(phi))

    G = np.array((x, y)).transpose()
    return G


def test_min_num(random_walk):
    # Set maximum allowable error
    epsilon = 20.0

    G_pp = polyprox.min_num(random_walk, epsilon)
    G_rdp = rdp.rdp(random_walk, epsilon)

    assert np.array_equal(G_rdp, G_pp)


def test_min_e(random_walk):
    # Number of points by which to approximate
    m = 10
    G_pp = polyprox.min_e(random_walk, m)
