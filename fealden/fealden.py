#!/usr/bin/env python3

import argparse
import multiprocessing
import re
import time
import timeit

from . import seed, sensor

# __all__ = ["Fealden", "Fold", "Node", "Seed", "Sensor", "Structure"]

BINDING_STATE = {"DS": 0, "SS": 1}
verbose = False

# Set seed graph patterns from literature
SEED_GRAPHS = {
    0: [
        "Seed Graph 1:",
        "2",
        "1 0 2",
        "2 1 3 3 5",
        "3 2 2",
        "5 2 0",
        "Seed Graph 2:",
        "2",
        "1 0 2",
        "2 1 3 3 0",
        "3 2 2",
        "Seed Graph 3:",
        "4",
        "1 0 2",
        "2 1 3 7 9",
        "3 2 4",
        "4 3 5 5 7",
        "5 4 4",
        "7 4 2",
        "9 2 0",
    ],
    1: [
        "Seed Graph 1:",
        "5",
        "1 0 2",
        "2 1 3 3 5",
        "3 2 2",
        "5 2 0",
        "Seed Graph 2:",
        "7",
        "2 1 3 11 0",
        "3 2 4",
        "4 3 5 5 7",
        "5 4 4",
        "7 4 6",
        "6 7 9 9 11",
        "9 6 6",
        "11 6 2",
        "Seed Graph 3:",
        "3",
        "1 0 2",
        "2 1 3 3 0",
        "3 2 2",
    ],
}


def __main__() -> None:
    """
    __main__() begins the program, parses, and validates the command line arguments.

    The program is run by typing "python ./fealden.py STRING INT' where:
        STRING is a string, in quotes, which consists of A, T, C, G, a, c, t, or g
        and which represents the recognition sequence
        INT is an integer, either 0 or 1, which denotes the binding state of the
        recognition sequence. 0 stands for DS, 1 stands for SS.

    Parameters:
        None
    Returns:
        Nothing
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "recSeq",
        type=str,
        help="The sequenced recognized by your target \
                represented as a string comprised of the letters 'a', 'A', 't', \
                'T', 'c', 'C', 'g', and 'G'.",
    )
    parser.add_argument(
        "bindingState",
        type=int,
        help="The state of the sequence when bound to the target.\
                \n This is 0 if your target binds to a double stranded sequence\
                and 1 if it binds to a single stranded sequence.",
    )
    parser.add_argument(
        "-ms",
        "-max_size",
        type=int,
        help="The maximum number of bases allowed in your sensor.\
                This number must be greater than 20 and should probobly be less \
                than 50.",
        default=50,
    )
    parser.add_argument(
        "-sps",
        "-sens_per_seed",
        type=int,
        help="The minimum number of potential sensors per seed graph \
                structure. \n This number should be increased if you\
                are not getting enough results. (Start by increasing it to 2000, \
                and increase from there if you are still unsatisfied.)",
        default=500,
    )
    parser.add_argument(
        "-i",
        "--interactive",
        action="store_true",
        help="Interactive Mode, no file output",
    )
    parser.add_argument(
        "-out",
        "-output_file",
        type=str,
        help="The output file to store results.\n Results are written in CSV format.",
        default=f'{time.strftime("%Y%m%d-%H%M%S")}-results.csv',
    )
    parser.add_argument(
        "-v",
        "-verbose",
        action="store_true",
        help="Output information when each thread starts and completes operation.",
    )
    # Up next: Binding affinity tuning
    # Up next: Anticipated target concentration tuning

    args = parser.parse_args()
    invalidChars = re.compile("[^atgc]", re.IGNORECASE)
    if invalidChars.search(args.recSeq):
        print(
            "Invalid recognition sequence.\n \
                You must enter a sequence consisisting only of a,c,t,g,A,C,T,G."
        )
        exit(0)
    if args.bindingState != 0 and args.bindingState != 1:
        print("Invalid binding state. Argument must be 0 or 1. See -h for help.")
        exit(0)
    if args.ms <= 20:
        print("Maximum sensor size is too low, it must be greater than 20.")
        exit(0)
    global verbose
    verbose = args.v
    Fealden(
        args.recSeq.lower(),
        args.bindingState,
        args.ms,
        args.sps,
        args.interactive,
        args.out,
    )


def generate_sensor(
    seed: seed.Seed, recSeq: str, numPossSen: int, core: int
) -> list[sensor.Sensor]:
    """
    generate_sensor() gens a # of possible sensors and returns a list of valid sensors.

    Parameters:
        seed        <-- an object of the 'Seed' class, the seed graph for the sensor
        recSeq      <-- a String, the recognition sequence
        numPossSen  <-- an integer, the number of possible sensors to be generated
        core        <-- an integer, the ID of the core in which this process is running

    Retuns:
        sensors     <-- list of objecfs of the class 'Sensor'
    """
    # global verbose
    # if verbose:
    #     print("Starting: %s, core %d" % (seed.name, core))

    sensors = []
    version = 0
    minScore = 0

    while version < numPossSen:
        version += 1
        sen = seed.build_sensor(core, version, recSeq)

        # only keep good sensors
        if sen is None:
            continue
        if sen.score >= minScore:
            sensors.append(sen)
    # print sorted(sensors)

    # if verbose:
    #     print("Completed: %s, core %d" % (seed.name, core))
    return sensors


# *************************************************************************************
# Generating a Fealden object auto-runs all non-interactive parts of the program.
# *************************************************************************************


class Fealden:

    """
    __init__() the constructor for Fealden objects. The time it takes to run the
    program is printed to the standard out.
    To run this method sequentially, uncomment the 8th line in the method body.

    Parameters:
        recSeq         <-- a string, the recognition sequence.
        bindingState   <-- an integer, 0 or 1, representing the binding state of the
                           recognition sequence.
        maxSensorSize  <-- an integer, the max number of bases the user would like
                           in their sensor.
        minSensPerSeed <-- an integer, the minimum number of potential sensors to be
                           generated per seed graph.
        interactive    <-- a bool for interactive mode to store results in output
                           attribute, rather than write to a csv file.
        outputfile     <-- a string, filename to store results in.
    Returns:
        an object of the class Fealden
    """

    def __init__(
        self,
        recSeq: str,
        bindingState: int,
        maxSensorSize: int,
        minSensPerSeed: int,
        interactive: bool,
        outputfile: str,
    ) -> None:
        """Initialize new Fealden instance."""
        self.recSeq = recSeq
        self.bindingState = bindingState
        self.maxSensorSize = maxSensorSize
        self.outputfile = outputfile

        # Pulled seed file constructs into program to reduce file reads and
        # remove file dependencies
        seeds = self.parse_seed_file(SEED_GRAPHS[int(bindingState)])

        # recommendedSensPerSeed = 10 * len(recSeq) *\
        #     ((50 - maxSensorSize) if (maxSensorSize < 40) else (10))
        posSensPerSeed = minSensPerSeed

        # FIXME: Implement recommendedSensPerSeed
        # recommendedSensPerSeed if \
        # recommendedSensPerSeed < minSensPerSeed else minSensPerSeed

        timeZero = timeit.default_timer()
        # '''
        numProcess = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(numProcess)
        seedSensPerProcess = posSensPerSeed / numProcess

        tasks = []
        sensors: dict[str, sensor.Sensor] = {}

        i = 0
        while i < numProcess:
            i += 1
            tasks.extend([(s, self.recSeq, seedSensPerProcess, i) for s in seeds])

        for t in tasks:
            pool.apply_async(
                generate_sensor,
                t,
                callback=lambda result: sensors.update(
                    [(r.seq, r) for r in result]  # type: ignore[attr-defined]
                ),
            )
        pool.close()
        pool.join()

        s = sorted(sensors.values(), key=lambda sen: sen.score)

        if not interactive:
            try:
                f = open(self.outputfile, "w")
            except OSError:
                print("Unable to open " + self.outputfile)
                self.outputfile = "recent-results.csv"
                f = open(self.outputfile, "w")

            f.write(sensor.Sensor.csv_header() + "\n")
            for sen in s:
                f.write(str(sen) + "\n")
            f.close()

            print("Stored " + str(len(s)) + " result(s) in " + self.outputfile)
            print("Took " + str(timeit.default_timer() - timeZero) + " seconds")
        else:
            output_list = []
            output_list.append(sensor.Sensor.csv_header())
            for sen in s:
                output_list.append(str(sen))
            self.output = output_list

    def parse_seed_file(self, lines: list[str]) -> list[seed.Seed]:
        """
        parse_seed_file() is a simple method for parsing the seedGraph file.
        A seed graph file looks like this:
            Seed Graph 1:   <- The following information is for the first seed graph
            2               <- The name of the node which contains the recognition seq
            2 1 3 3 5
            3 2 2
            5 2 4
            4 5 7 7 0
            7 4 4
            Seed Graph 2:
            4
            1 0 2
            2 1 3 3 5
            3 2 2
            5 2 4
            4 5 7 7 0
            7 4 4
            .
            .
            .
        See the comment for the "make_graph" method of the seed class for an explenation
        of how the data given corrilates to a graph.

        This method breaks the data into chunks, one per seed graph, where each chunk is
        a list containing all the info about that particular graph. Then each list
        is used to generate a new seed graph using the class Seed. A list of Seed objs
        is returned.

        Parameters:
            lines <-- a list of str, each str is a line in the seed graph data file

        Reutrns:
            seeds <-- a list of objects of the class "Seed"
        """

        seeds = []
        i = 0
        li = lines[i]
        seedNum = 0
        while i < (len(lines)):
            i += 1
            li = lines[i]
            graphData = []
            recNodeName = "-1"
            seedNum += 1
            while i < (len(lines)) and lines[i].split()[0] != "Seed":
                li = lines[i].strip()
                if len(li) == 1:
                    recNodeName = li
                    i += 1
                    continue
                graphData.append(li)
                i += 1
            seeds.append(
                seed.Seed(
                    graphData,
                    recNodeName,
                    self.recSeq,
                    self.bindingState,
                    str("Graph " + str(seedNum)),
                    self.maxSensorSize,
                )
            )
        return seeds


# run the program from the command line
if __name__ == "__main__":
    __main__()
