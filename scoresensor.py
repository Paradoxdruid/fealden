#!/usr/bin/env python3

import argparse, re, os, subprocess, shutil, tempfile
import seed, sensor


bindingstate = 1
tempdir = 'tmp'


"""
Main function
Parses arguments and prints sensor info to stdout in CSV format
"""
def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('senSeq', type = str,
                        help = 'The full sensor sequence.')
    parser.add_argument('bindLoc', type = int,
                        help = 'Index of the binding sequence within the sensor sequence, beginning from 0.')
    parser.add_argument('bindLen', type = int,
                        help = 'Length of the binding sequence within the sensor sequence.')

    args = parser.parse_args()
    invalid = re.compile('[^atgc]', re.IGNORECASE)
    if invalid.search(args.senSeq):
        print('Invalid sensor sequence.\nYou must enter a sequence consisting only of a, t, g, c, A, T, G, and C.')

    try:
        sen = build_custom_sensor(args.senSeq, args.bindLoc, args.bindLen)
        if sen.score >= 0:
            print(sen)
        else:
            print('Not a functional sensor: code ' + str(sen.score))
    except:
        print('An error has occurred')
        raise


"""
Build a sensor for evaluation
params: sequence string of sensor
        binding location, integer index of the binding region within the sensor
        binding length, integer length of binding region
return: a scored sensor
"""
def build_custom_sensor(sequence_string, binding_location, binding_length):
    sequence_list = list(sequence_string.upper())
    seq_file_name = '-'.join([sequence_string[:4], str(binding_location), str(binding_length), 'eval'])

    try:
        seqfile = open(os.path.join(tempdir, seq_file_name), 'w')
        seqfile.write(str(sequence_list))
        seqfile.close()
    except:
        print('Cannot create sequence file')
        raise

    try:
        subprocess.check_call(seed.command + [seq_file_name], cwd=tempdir, stdout=open("/dev/null"))
    except:
        print('Call to ' + seed.command + ' failed')
        raise

    recSeq = {'start': binding_location + 1,
              'end':   binding_location + 1 + binding_length}
    respSeq = {'start': -1, 'end': -1}

    foldfile = open(os.path.join(tempdir, seq_file_name + '.ct'), 'r')
    sen = sensor.Sensor(foldfile, recSeq, respSeq, bindingstate, 'eval', sequence_string[binding_location:binding_location + binding_length])
    foldfile.close()
    return sen


# Entry point of script
if __name__ == '__main__':
    tempdir = tempfile.mkdtemp()
    __main__()
    shutil.rmtree(tempdir)

