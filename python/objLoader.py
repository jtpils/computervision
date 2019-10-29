#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

class ObjLoader(object):
    def __init__(self, file_name):
        self.vertices = []

        try:
            f = open(file_name)
            for line in f:
                if line[:2] == 'v ':
                    index1 = line.find(' ') + 1
                    index2 = line.find(' ', index1 + 1)
                    index3 = line.find(' ', index2 + 1)

                    vertex = (float(line[index1:index2]),
                              float(line[index2:index3]),
                              float(line[index3:-1]))

                    vertex = [round(vertex[0], 2),
                              round(vertex[1], 2),
                              round(vertex[2], 2)]

                    self.vertices.append(vertex)

                elif line[0] == 'f':
                    pass

            f.close()

            self.vertices = np.array(self.vertices)
        except IOError:
            print('.obj file not found.')
