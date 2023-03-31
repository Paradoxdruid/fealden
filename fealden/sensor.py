from . import fold


class Sensor:

    """
    ---------------------------------------------------------------
    Sensor is a structure to hold and interpret the results from
    a unafold query, along with some other information. A Sensor
    object has a sequence, a recognition sequence and
    its location in the overall sequence, a list of Folds,
    a tagging location, the desired state of the recognition
    sequence, the name of the seed graph that gave rise to
    this sensor, and an overall score.
    -----------------------------------------------------------------
    """

    def __init__(
        self,
        data_file: tuple[str, list[dict[str, float | list[list[int]]]]],
        rec_seq: dict[str, int],
        resp_seq: dict[str, int],
        des_rec_seq_state: int,
        seed_name: str,
        base_seq: str,
    ):
        """
        This is the constructor for Sensor.
        Paramaters:
                data_file <- a File object (this should be the .ct
                                            returned from a unafold query)
                rec_seq   <- a dict of the form {'start': n, 'end': p}, where n,p are
                            integers and represent, respectivly, the starting and
                            ending locaiton of the recognition sequence within the
                            overall sensor's sequence.
                resp_seq <- a dict of the form {'start': n, 'end': p}, where n,p are
                            integers and represent, respectivly, the starting and
                            ending locaiton of the recognition response sequence
                            within the overall sensor's sequence. n and p are both -1
                            if the recognition sequence is single stranded.
                des_rec_seq_state <- An integer, either 0 or 1, representing the state
                            in which the recognition sequence binds to the target.
                            0 represents double stranded, 1 represents single stranded.
                seed_name <- An integer, this is a simple tag to represent which
                            graph gave rise to this sensor.
        """
        self.seed_name = seed_name
        self.rec_seq = rec_seq
        self.resp_seq = resp_seq
        self.des_rec_seq_state = des_rec_seq_state

        # Data file is passed to interpret_data, actually list
        (self.seq, self.folds) = self.interpret_data(data_file)
        self.on_conc = 0
        self.off_conc = 0
        self.noise_conc = 0
        self.on_to_off_dist = 0
        self.base_seq = base_seq
        (self.tag_loc, self.score) = self.get_tag_and_score()

    def interpret_data(
        self, data: tuple[str, list[dict[str, float | list[list[int]]]]]
    ) -> tuple[str, list[fold.Fold]]:
        """
        interpret_data takes data from a the .ct file which has
        been parsed into a list of lines. It returns a tuple.
        The first value is the sequence, represented as a string
        of lowercase letters. The second value is a list of Folds.

        Parameters:
            data   <-- A list of lines from the .ct file output by unafold.
        Returns:
            (seq, folds) <-- the tuple described above.

        """
        seq, structure_data = data[0], data[1]
        # (seq, structureData) = self.simplify_input(lines)
        folds: list[fold.Fold] = []

        for each in structure_data:
            this_fold = fold.Fold(
                each["bps"], each["deltaG"], self.rec_seq  # type: ignore
            )
            folds.append(this_fold)

        return (seq, folds)

    def simplify_input(
        self, lines: list[str]
    ) -> tuple[str, list[dict[str, float | list[int]]]]:
        """
        simplify_input takes the .ct file, represented as as list
        of lines, and distills from it the information we care
        about into a useable form. The information returned is
        a tuple. The first value is the sequence represented as
        a string. The second value is a list of dictionaries.
        The list consists of all the folds in the data file. The
        dictionary contains all the relevent information we have
        about each fold up to this point. See the following example
        for clarification.

        ('attcgtgcatggtcaatcttacgttacgacggcccattcaaa' ,
         [{'deltaG': -13.4 ,        <- The delta G for the first fold in the list
          'bps'   : [[1, 0],    <- base pair 1 is bound to nothing
                     [2, 0],
                     [3, 22],   <- Base pair 3 is bound to base pair 22
                     .
                     .
                     .
                     [44, 30]]
          }
          {'deltaG': -12.3 ,    <- The second fold in the list has a delta G of -12.3
           'bps'   : [[1, 23],  <- base pair 1 is bound to base pair 23
                     .
                     .
                     .
                     [44, 0]]   <- Base pair 44 is not bound
          }
          .
          .
          .
         ]
        )


        Parameters:
            lines <-- a list of lines from the .ct file
        Returns:
            (sequence, stuructureData) -- the tuple described above
        """
        structure_data = []
        fold_index = -1
        fold_size = int(lines[0].split()[0])
        sequence = []

        for i, v in enumerate(lines):
            if i % (fold_size + 1) == 0:
                # This line holds a deltaG value for a new structure
                deltaG = float(v.split()[3])
                structure_data.append({"deltaG": deltaG, "bps": []})
                fold_index += 1
            else:
                # This line holds information about the structure we're
                # currently working on
                temp = v.split()
                strlist = temp[0:1] + temp[4:5]
                structure_data[fold_index]["bps"].append(  # type: ignore
                    [int(x) for x in strlist]
                )
                # we've just appended the base pair number, and the number of the
                # base pair it is bound to
                if i < fold_size + 1:
                    sequence.append(temp[1])
        return ("".join(sequence).lower(), structure_data)  # type: ignore

    def get_tag_and_score(self) -> tuple[int, float]:
        """
        get_tag_and_score returns a number for the sensor, based on an arbitrary fitness
        scale, this serves as a proxy for our estimation of how well the sensor
        will work. The greater the number, the higher the estimation of the
        probability for success. When the value -1 is returned, the sensor was
        determined to be invalid entirely.

        Parameters:
            None
        Returns:
            A tuple, (x,y). x is an integer, the location on the sensor's sequence upon
            upon which the tag should be placed. y is a floating point number, the score
            of the sensor.
        """

        DELTA_G_MAX_DIFFERENCE = 5
        if len(self.folds) <= 1:
            # 'Only one fold'
            return (0, -1)
        if self.folds[1].deltaG - DELTA_G_MAX_DIFFERENCE > self.folds[0].deltaG:
            # "First two folds have delta Gs which are too disparate."
            return (0, -2)
        if (
            len(self.folds) > 2
            and self.folds[2].deltaG - DELTA_G_MAX_DIFFERENCE > self.folds[1].deltaG
        ):
            # "Delta Gs of 2nd and 3rd folds are disparate."
            if self.folds[0].rec_seq_state == self.folds[1].rec_seq_state:
                # "Recognition sequence is in the same state in the first
                # two folds."
                return (0, -3)
            if (
                self.folds[0].rec_seq_state != self.des_rec_seq_state
                and self.folds[1].rec_seq_state != self.des_rec_seq_state
            ):
                # "In neither of the first two folds is the recognition
                # sequence in the desired state."
                return (0, -4)
        if self.folds[0].deltaG > -2 or self.folds[0].deltaG < -50:
            # "The first has a delta G which is out of range."
            return (0, -5)
        # sensor has passed triage criteria
        # compute validity based on criteria requiring distance
        score_data = self.get_tagging_information()

        if score_data == 0:
            return (0, -6)
        (
            self.tag_loc,
            self.on_conc,
            self.off_conc,
            self.noise_conc,
            self.wrong_conc,
            self.fuzzy_conc,
            self.on_to_off_dist,
        ) = score_data  # type: ignore
        # if sensor is valid, get optimal tagging scenario
        # score sensor based on optimal tagging scenario
        return (self.tag_loc, self.calculate_score())

    def calculate_score(self) -> float:
        """Calculate a score for the sensor based on optimal desired scenario.

        This function looks to maximize change in calculated on/off distance
        as well as minimize difference in calculated concentration of ON and OFF state.
        """

        # Improved scoring function:
        # First term ranges from 0.5 (good) to 1 (bad), with good being Ks near 1
        # abs(2A-B) / (A+B) when A == B is 0.5, when B >> A tends to 1
        #
        # Second term, for current distance algorithm, tends 0.5 (good) to 1 (bad),
        # with larger distance favored; observed distances range from 10 to ~ 25
        # 1/10 (bad) = 0.1, 1/50 (good) = 0.04
        # Multiplied by 10 for rough parity in weighting, this could be adjusted

        SIGNAL_GAIN_WEIGHT = 10

        return (
            abs(2 * self.on_conc - self.off_conc) / (self.on_conc + self.off_conc)
        ) + (SIGNAL_GAIN_WEIGHT * (1 / self.on_to_off_dist))

    def get_onstate_and_wrong(
        self, MAX_ON_DIST: int, distances: list[int]
    ) -> tuple[list[tuple[int, float]], float]:
        """Helper for finding tagging information, finding onState and concWrong."""
        on_state_info = []
        conc_wrong: float = 0.0
        for i, d in enumerate(distances):
            # the fold to which this distance is refering
            curr_fold = self.folds[i]
            # distance is less than or equal to on dist (ie this is a on
            # fold)
            if MAX_ON_DIST - d >= 0:
                if self.des_rec_seq_state == curr_fold.rec_seq_state:
                    # This is a on position and the sensor will bind the
                    # target
                    on_state_info.append((d, curr_fold.conc))
                else:  # this is sending the opposite of the desired signal
                    conc_wrong += curr_fold.conc
                    # this is only to speed it up when we don't want any
                    # noise
                    break

        return on_state_info, conc_wrong

    def get_offstate_wrong_and_fuzzy(
        self,
        MAX_ON_DIST: int,
        MIN_OFF_CHANGE: int,
        distances: list[int],
        weighted_avg_on_dist: float,
    ) -> tuple[list[tuple[int, float]], float, float]:
        """Helper for finding tagging information, finding offState,
        concFuzzy and concWrong."""
        off_state_info = []
        conc_wrong = 0.0
        conc_fuzzy = 0.0
        for i, d in enumerate(distances):
            curr_fold = self.folds[i]
            if MAX_ON_DIST - d >= 0:
                continue  # we've already delt with these folds
            if MIN_OFF_CHANGE <= d - weighted_avg_on_dist:
                if self.des_rec_seq_state != curr_fold.rec_seq_state:
                    # This is a off position and the sensor will not bind
                    # the target
                    off_state_info.append((d, curr_fold.conc))
                else:  # this is sending the wrong signal
                    conc_wrong += curr_fold.conc
                    # this is only to speed it up when we don't want any
                    # noise
                    break
            # the tag distance is not close enough nor far enough to be off
            else:
                conc_fuzzy += curr_fold.conc
                # this is only to speed it up when we don't want any noise
                break
        return off_state_info, conc_wrong, conc_fuzzy

    def get_tag_locations(
        self, MAX_ON_DIST: int, MIN_OFF_CHANGE: int
    ) -> list[tuple[int, list[int]]]:
        tag_locs = []
        # get potential tagging locations and their distances in all the
        # various folds
        for i, v in enumerate(self.seq):
            if (
                v.lower() == "t"
                and (i + 1 < self.rec_seq["start"] or i + 1 > self.rec_seq["end"])
                and (i + 1 < self.resp_seq["start"] or i + 1 > self.resp_seq["end"])
            ):
                distances = [f.get_distance(1, i + 1) for f in self.folds]
                smallest_dist = min(distances)
                if (
                    smallest_dist <= MAX_ON_DIST
                    and max(distances) - smallest_dist >= MIN_OFF_CHANGE
                ):
                    tag_locs.append((i + 1, distances))
        return tag_locs

    def get_tagging_information(
        self,
    ) -> int | tuple[int, float, float, float, float, float, float]:
        """
        get_tagging_information() finds the optimal tagging situation, and returns some
        information with which it is associated. The information returned is labeled in
        the "Returns:" section below and should be interpretd as follows:
        position               -- the location, expressed as an integer, of the base
                                    to be tagged
        onConc                 -- the (relative) concentration of sensor that's 'on.'
        offConc                -- the (relative) concentration of sensor that's 'off.'
        noiseConc              -- the (relative) concentration of sensors that are
                                    'wrong' or 'fuzzy.'
        concWrong              -- the relative concentration of sensor that is on/off
                                    when it should be the opposite. (ie it's in the
                                    wrong state.)
        concFuzzy              -- the relative concentration of sensor in which the
                                    tags are neither close enough to be 'on' nor far
                                    enough away from the 'on' states to be truly 'off.'
        weightedAvgOnToOffDist -- the average distance of the on states
                                    less the average distance of the off states.
                                    All averages are weighted based on the
                                    concentrations of the various on and off states.


        Parameters:
            None

        Returns:
            (position, onConc, offConc, noiseConc,
            concWrong, concFuzzy, weightedAvgOnToOffDist)
        """
        MAX_ON_DIST = 12
        MIN_OFF_CHANGE = 10
        tag_locs = self.get_tag_locations(MAX_ON_DIST, MIN_OFF_CHANGE)

        max_avg_delta_on_to_off = 0
        # determine if this would make a good sensor if tagged in each possible
        # location
        for t in tag_locs:
            on_conc = 0.0
            off_conc = 0.0
            noise_conc = 0.0
            conc_wrong = 0.0
            conc_fuzzy = 0.0

            # position == the loation of the tag
            # distances == the physical distance between the tag
            #             and the first position in each fold
            (position, distances) = t

            on_state_info, new_conc_wrong = self.get_onstate_and_wrong(
                MAX_ON_DIST, distances
            )
            conc_wrong += new_conc_wrong

            if on_state_info == []:  # no on states
                continue

            # total concentration of all the on states
            on_conc = sum([j for (i, j) in on_state_info])

            weighted_avg_on_dist = sum([i * (j / on_conc) for (i, j) in on_state_info])

            (
                off_state_info,
                more_concWrong,
                new_concFuzzy,
            ) = self.get_offstate_wrong_and_fuzzy(
                MAX_ON_DIST, MIN_OFF_CHANGE, distances, weighted_avg_on_dist
            )
            conc_wrong += more_concWrong
            conc_fuzzy += new_concFuzzy

            noise_conc = conc_fuzzy + conc_wrong

            if off_state_info == []:  # no off states
                continue

            # the concentration of all the off states
            off_conc = sum([j for (i, j) in off_state_info])

            # if noiseConc*10 > offConc + onConc or\
            #    concWrong*10 > offConc or\
            #    concWrong*10 > onConc or\

            if noise_conc > 0 or off_conc * 10 < on_conc or on_conc * 10 < off_conc:
                # "too much noise, too many are wrong, or ratios are off"
                continue
            #
            weighted_avg_on_to_off_dist = sum(
                [
                    (i - weighted_avg_on_dist) * (j / off_conc)
                    for (i, j) in off_state_info
                ]
            )

            # CHANGE to check for best overall score
            if max_avg_delta_on_to_off < weighted_avg_on_to_off_dist:
                return (
                    position,
                    on_conc,
                    off_conc,
                    noise_conc,
                    conc_wrong,
                    conc_fuzzy,
                    weighted_avg_on_to_off_dist,
                )

        return 0

    @staticmethod
    def csv_header() -> str:
        """
        csv_header() generates the header to a CSV (comma separated values) file that
        describes the values returned by csv_line().

        Parameters:
            None
        Returns:
            A string, the comma delimited header for a CSV file containing information
            about sensors.
        """
        return ",".join(
            [
                "Sequence",
                "Score",
                "Seed Name",
                "Tag Location",
                "Conc On",
                "Conc Off",
                "Conc Off to On Ratio",
                "Conc Noise",
                "Conc Wrong",
                "Conc Fuzzy",
                "On to Off Dist (nm)",
                "Length",
                "Num Folds",
                "Original Sequence",
            ]
        )

    def csv_line(self) -> str:
        """
        csv_line() generates a string representing this sensor, which can be inserted
        into a CSV (comma separated values) file.

        Parameters:
            None
        Returns:
            A string, the string representation of a sensor, comma delimited.
        """
        return ",".join(
            [
                self.seq,
                str(self.score),
                str(self.seed_name),
                str(self.tag_loc),
                str(self.on_conc),
                str(self.off_conc),
                str(self.off_conc / self.on_conc),
                str(self.noise_conc),
                str(self.wrong_conc),
                str(self.fuzzy_conc),
                str(self.on_to_off_dist),
                str(len(self.seq)),
                str(len(self.folds)),
                self.base_seq,
            ]
        )

    def __repr__(self) -> str:
        """
        __repr__() generates the string representation of a sensor. It is all
        the information one might want to know about a sensor.


        Parameters:
            None
        Returns:
            A string, the string representation of a sensor.
        """

        return self.csv_line()
