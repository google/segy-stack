import pytest

from segystack import Grid, UTMZone

class TestGrid:
    def test_init(self):
        g = Grid()
        assert(g.inline_min == 0)
        g.inline_min = 12
        assert(g.inline_min == 12)
        assert(g.inline_max == 0)
        g.inline_max = 120
        assert(g.inline_max == 120)
        print(g)

    def test_types(self):
        c = Grid.Cell()
        assert(c.x_coordinate == 0.0)
        c.x_coordinate = 3.1
        assert(c.x_coordinate == 3.1)
        print(c)

        utm = UTMZone()
        utm.letter = "H"
        utm.number = 38
        print(utm)

if __name__ == "__main__":
    t = TestGrid()
    t.test_init()
    t.test_types()
