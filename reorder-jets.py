#! /usr/bin/python

# Experimental!
# Totally unsuitable for any purpose other than causing you pain!


# Run this as
#  python reorder-intakes.py foo.craft
# It will output 
#   foo.craft-reordered 
# with intakes and engines reordered to try to make for more symmetric
# flameout.  Check it out, rename it to foo2.craft or whatever, and 
# load it in the game.  Be very careful, this isn't robust in any way.

from __future__ import division

import sys
import re
import jets

class part(object):
    def __init__(self, stream):
        """
        Given a stream that just read "PART", read everything
        from the next bracket to the final bracket and keep track
        of whether it's an intake or air-breathing engine.

        Assumption: not much happens in one line.
        """
        brackets = 0
        self.lines = []

        self.intakeType = None
        self.engineType = None

        for line in stream:
          self.lines.append(line)
          if re.search('}', line):
              brackets -= 1
              if brackets == 0:
                # this part is done!
                return
          elif re.search('{', line):
              brackets += 1
          else:
            for x in jets.engines:
              if re.search('part = %s[0-9]*' % x.partName, line):
                  self.engineType = x
            for x in jets.intakes:
              if re.search('part = %s[0-9]*' % x.partName, line):
                  self.intakeType = x

    def toStream(self, stream):
        stream.write("PART\n")
        stream.writelines(self.lines)

class ship(object):
    @staticmethod
    def reorderFile(name):
        with open(name) as f:
            vessel = ship(f)
        vessel.summarize(sys.stdout)
        with open(name + "-reordered", 'w') as f:
            f.writelines(vessel.prefix)
            for part in vessel.reorder():
              part.toStream(f)

    def _append(self, part, type, dict):
        if type is None:
          return False
        if type not in dict:
          dict[type] = [part]
        else:
          dict[type].append(part)
        return True

    def __init__(self, stream):
        # map type => list of part
        self.engines = dict()
        self.intakes = dict()
        # list of other parts
        self.otherParts = []
        # list of all parts in original order
        self.allParts = []
        # list of lines that come before the first PART
        self.prefix = []
        for lineno, line in enumerate(stream):
            if not re.search('PART', line):
              if len(self.allParts) == 0:
                self.prefix.append(line)
              else:
                print ("uninterpreted text at line %d: %s" % (lineno, line))
            else: # line == 'PART':
              nextPart = part(stream)
              self.allParts.append(nextPart)
              if self._append(nextPart, nextPart.engineType, self.engines):
                  continue
              if self._append(nextPart, nextPart.intakeType, self.intakes):
                  continue
              self.otherParts.append(nextPart)

    def summarize(self, stream):
        # for debugging, print the characteristics
        stream.write("%d parts\n" % len(self.allParts))
        stream.write("%d engine types\n" % len(self.engines))
        for eType, eList in self.engines.iteritems():
            stream.write("  %d %s\n" % (len(eList), eType.name))
        stream.write("%d intakes\n" % len(self.intakes))
        for iType, iList in self.intakes.iteritems():
            stream.write("  %d %s\n" % (len(iList), iType.name))


    def reorder(self):
        # dump in all the other parts first.
        parts = [ x for x in self.otherParts ]

        # Iterate to get the last key/value pair, or the turbojets.
        for jettype in self.engines:
          if jettype == jets.turbojet:
            break
        jetlist = self.engines[jettype]

        # for each engine, set up a list of parts
        numEngines = len(jetlist)
        partsByEngine = [ [x] for x in jetlist ]
        for iType, iList in self.intakes.iteritems():
          numIntakes = len(iList)
          for i, intake in enumerate(iList):
              # allocate intakes in round-robin
              partsByEngine[i % numEngines].append(intake)

        # add each engine and its intakes to the parts list.
        # make sure to reverse the order, so the intakes go before the engine
        for eParts in partsByEngine:
          parts.extend(reversed(eParts))

        # add all the other engines to the parts list (if any), ordered by type
        for t in self.engines:
          if t != jettype:
            parts.extend(self.engines[t])

        return parts

for name in sys.argv[1:]:
    ship.reorderFile(name)
