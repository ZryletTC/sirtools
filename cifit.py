import numpy as np
import argparse


def read_mch_params(lines, n_expected, specify_vary=False):
    """
    Reads a line from the given iterator `lines` and returns the parameters.

    If `specify_vary` is True, it also reads whether the parameters should be
    varied from the next line. Otherwise, it assumes the user will be asked
    about constraints later.
    """

    param_guess = list(map(float, next(lines).split()))
    assert len(param_guess) == n_expected

    if specify_vary:
        param_vary = list(map(bool, next(lines).split()))
        assert len(param_vary) == n_expected
    else:
        param_vary = None

    return param_guess, param_vary


def read_mch_file(filename, specify_vary=False):
    print(f"Specifying vary: {specify_vary}")
    with open(filename, 'r') as file:
        lines = file.readlines()

    lines = iter([line.strip() for line in lines
                  if line.strip() and not line.startswith('#')])

    # Title line
    mech_title = next(lines)
    print(f"Mechanism title: {mech_title}")

    # Read number of sites and processes
    n_sites, n_procs = map(int, next(lines).split())

    # Read T1 rate(r1), M0 and Minf guesses, and check against number of sites
    r1_guess, r1_vary = read_mch_params(lines, n_sites,
                                        specify_vary=specify_vary)
    minf_guess, minf_vary = read_mch_params(lines, n_sites,
                                            specify_vary=specify_vary)
    m0_guess, m0_vary = read_mch_params(lines, n_sites,
                                        specify_vary=specify_vary)

    # Read parameters for each process
    processes = []
    for i in range(n_procs):
        # First line of each process block specifies the rate constant k
        # and whether it should be varied or not
        if specify_vary:
            k_guess_str, k_vary_str = next(lines).split()
            k_guess = float(k_guess_str)
            k_vary = bool(k_vary_str)
        else:
            k_guess = float(next(lines))
            k_vary = None

        # Initialize diagonal matrix for the process
        m = np.diag(np.ones(n_sites))

        # Read the nonzero off-diagonals associated with the process,
        # and add them to the matrix
        n_off_diagonals = int(next(lines))
        for j in range(n_off_diagonals):
            row_str, col_str, val_str = next(lines).split()
            row = int(row_str)  # Row of off-diagonal element
            col = int(col_str)  # Column of off-diagonal element
            val = -float(val_str)  # Invert sign as per original code
            m[row, col] = val

        process = {
            'k_guess': k_guess,
            'k_vary': k_vary,
            'matrix': m
        }
        processes.append(process)

    # If not specifying vary, ask the user for each parameter
    if not specify_vary:
        def ask_vary(param_name, guesses):
            print(f"Should {param_name} parameters vary during fitting?")
            vary = []
            for i, guess in enumerate(guesses):
                while True:
                    prompt = f"  Site {i+1} ({param_name} = {guess}): [y/n] "
                    resp = input(prompt).strip().lower()
                    if resp in ('y', 'n', 'yes', 'no'):
                        vary.append(resp[0] == 'y')
                        break
                    else:
                        print("Please enter 'y' or 'n'.")
            return vary

        r1_vary = ask_vary("R1", r1_guess)
        minf_vary = ask_vary("Minf", minf_guess)
        m0_vary = ask_vary("M0", m0_guess)

        k_vary = ask_vary("k", [proc['k_guess'] for proc in processes])
        for i, proc in enumerate(processes):
            proc['k_vary'] = k_vary[i]

    mch_dict = {
        'title': mech_title,
        'n_sites': n_sites,
        'n_procs': n_procs,
        'r1_guess': r1_guess,
        'r1_vary': r1_vary,
        'minf_guess': minf_guess,
        'minf_vary': minf_vary,
        'm0_guess': m0_guess,
        'm0_vary': m0_vary,
        'processes': processes
    }

    assert all(type(vary) is bool for vary in r1_vary + minf_vary + m0_vary), \
        "All vary parameters should be boolean values."

    return mch_dict


def read_data_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    lines = iter([line.strip() for line in lines
                  if line.strip() and not line.startswith('#')])

    # Read title line
    data_title = next(lines)
    print(f"Data title: {data_title}")

    # Read number of time points
    n_points = int(next(lines))

    # Initialize arrays for time points and data values
    time_points = []
    data_values = []

    for i in range(n_points):
        line = next(lines).split()
        # First element is the time point
        time_points.append(float(line[0]))
        # Remaining elements are the magnetization values
        data_values.append(list(map(float, line[1:])))

    # Transpose to get magnetizations per site
    magnetizations = np.array(data_values).T

    return {
        'title': data_title,
        'n_points': n_points,
        'time_points': time_points,
        'magnetizations': magnetizations
    }


def print_parameters(par_dict):
    print(f"Number of sites: {par_dict['n_sites']}")
    for i in range(par_dict['n_sites']):
        print(f"Site {i+1}:")
        print(f"  R1 relaxation rate: {par_dict['r1_guess'][i]:.4f}, "
              "vary: {par_dict['r1_vary'][i]}")
        print(f"  M_inf: {par_dict['minf_guess'][i]}, "
              "vary: {par_dict['minf_vary'][i]}")
        print(f"  M_0: {par_dict['m0_guess'][i]}, "
              "vary: {par_dict['m0_vary'][i]}")
        print()

    print(f"Number of processes: {par_dict['n_procs']}")
    for i, proc in enumerate(par_dict['processes']):
        print(f"Process {i+1}: k = {proc['k_guess']}, vary = {proc['k_vary']}")
        print("Process matrix:")
        print(proc['matrix'])
        print()


def parse_args():
    parser = argparse.ArgumentParser(description="CIFIT: A program for fitting"
                                     " selective inversion experiment data.")
    parser.add_argument("filename", help="Name of .dat and .mch files to"
                        " process, and of .out file to write.")
    parser.add_argument("--specify-vary", action="store_true",
                        help="Specify whether parameters should vary"
                             " during fitting inside the .mch file.")
    return parser.parse_args()


def main():
    args = parse_args()

    data_filename = f"{args.filename}.dat"
    mech_filename = f"{args.filename}.mch"
    out_filename = f"{args.filename}.out"

    # Check if files exist and are readable/writable
    try:
        with open(data_filename, 'r'), open(mech_filename, 'r'):
            pass
        print(f"Mechanism file: {mech_filename}")
        print(f"Data file: {data_filename}")
    except Exception as e:
        print(f"Error opening files: {e}")
        return 1

    try:
        with open(out_filename, 'w'):
            pass
        print(f"Output file: {out_filename}")
    except Exception as e:
        print(f"Error creating output file: {e}")
        return 1

    # Read mechanism file
    mch_dict = read_mch_file(mech_filename, specify_vary=args.specify_vary)
    print(mch_dict)

    # Read data file
    data_dict = read_data_file(data_filename)

    print_parameters(mch_dict)


if __name__ == "__main__":
    main()
