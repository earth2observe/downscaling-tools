#!/usr/bin/python

# e2o_dstools is Free software, see below:
#
# Copyright (c) Deltares 2005-2014
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.



# $Rev:: 904           $:  Revision of last commit
# $Author:: schelle    $:  Author of last commit
# $Date:: 2014-01-13 1#$:  Date of last commit

"""
e2o_resample_stack -- resample a stack of pcraster maps

Usage::

    -D dirname
    -O output dir name
    -C clone map (map to interpolate to)
    -M maxcpu
       maximum number of cpu's/cores to use (default = 4)

The script uses the pcraster resample program to resample the maps.
"""





import os, sys, shlex, time
import os.path
import glob
import getopt
import subprocess



def usage(*args):
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)



def runCommands(commands, maxCpu):
    """
    Runs a list of processes dividing
    over maxCpu number of cores.
    """

    def removeFinishedProcesses(processes):
        """ given a list of (commandString, process),
            remove those that have completed and return the result
        """
        newProcs = []
        for pollCmd, pollProc in processes:
            retCode = pollProc.poll()
            if retCode==None:
                # still running
                newProcs.append((pollCmd, pollProc))
            elif retCode!=0:
                # failed
                raise Exception("Command %s failed" % pollCmd)
            else:
                print "Command %s completed successfully" % pollCmd
        return newProcs

    processes = []
    for command in commands:
        command = command.replace('\\','/') # otherwise shlex.split removes all path separators
        proc =  subprocess.Popen(shlex.split(command))
        procTuple = (command, proc)
        processes.append(procTuple)
        while len(processes) >= maxCpu:
            time.sleep(.2)
            processes = removeFinishedProcesses(processes)

    # wait for all processes
    while len(processes)>0:
        time.sleep(0.5)
        processes = removeFinishedProcesses(processes)
    print "All processes in que (" + str(len(commands)) + ") completed."


def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'D:O:C:M:Nf')
    except getopt.error, msg:
        usage(msg)

    clone = "clone.map"
    Verbose=1
    inputdir = "inmaps"
    outputdir = "outmaps"
    mapstack = []
    maxcpu = 4
    force = False
    dirs = []

    for o, a in opts:
        if o == '-C': clone = a
        if o == '-D': inputdir = a
        if o == '-O': outputdir = a
        if o == '-f': force = False
        if o == '-M': maxcpu = int(a)


    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)



    allcmd = []
    for mfile in glob.glob(inputdir + '/*.[0-9][0-9][0-9]'):
        mstr = "resample --clone " + str(clone) + ' ' + mfile + " " + mfile.replace(inputdir,outputdir)
        if not os.path.exists(mfile.replace(inputdir,outputdir)):
            allcmd.append(mstr)
        else:
            print "skipping existing file: " + mfile.replace(inputdir,outputdir)
    runCommands(allcmd,maxcpu)




if __name__ == "__main__":
    main()
