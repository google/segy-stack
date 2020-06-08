#!/usr/bin/python
#
# Copyright 2020 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pytest

from segystack import StackFile, UTMZone

class TestGrid:
    def test_init(self):
        g = StackFile.Grid()
        assert(g.inline_min == 0)
        g.inline_min = 12
        assert(g.inline_min == 12)
        assert(g.inline_max == 0)
        g.inline_max = 120
        assert(g.inline_max == 120)
        print(g)

    def test_types(self):
        utm = UTMZone()
        utm.letter = "H"
        utm.number = 38
        print(utm)

if __name__ == "__main__":
    t = TestGrid()
    t.test_init()
    t.test_types()
