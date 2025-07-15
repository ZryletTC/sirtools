#!/usr/bin/env python3

import argparse

import numpy as np
from lmfit import Parameters, minimize, report_fit

VERSION = "0.1"  # TODO: Update to dynamic versioning from package


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
        param_vary = list(map(bool, map(int, next(lines).split())))
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
            k_vary = bool(int(k_vary_str))
        else:
            k_guess = float(next(lines))
            k_vary = None

        # Initialize diagonal matrix for the process
        # exch_mat = np.diag(np.ones(n_sites))
        exch_mat = np.zeros((n_sites, n_sites))

        # Read the nonzero off-diagonals associated with the process,
        # and add them to the matrix
        n_off_diagonals = int(next(lines))
        for j in range(n_off_diagonals):
            row_str, col_str, val_str = next(lines).split()
            row = int(row_str)  # Row of off-diagonal element
            col = int(col_str)  # Column of off-diagonal element
            val = -float(val_str)  # Invert sign as per original code
            exch_mat[row, col] = val

        process = {
            'k_guess': k_guess,
            'k_vary': k_vary,
            'matrix': exch_mat
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

    const_dict = {
        'title': mech_title,
        'n_sites': n_sites,
        'n_procs': n_procs,
        'matrices': []
    }

    # TODO: Remove 'guess' from key names
    pars_dict = {
        'r1_guess': r1_guess,
        'r1_vary': r1_vary,
        'minf_guess': minf_guess,
        'minf_vary': minf_vary,
        'm0_guess': m0_guess,
        'm0_vary': m0_vary,
        'k_guess': [],
        'k_vary': []
    }
    for proc in processes:
        pars_dict['k_guess'].append(proc['k_guess'])
        pars_dict['k_vary'].append(proc['k_vary'])
        const_dict['matrices'].append(proc['matrix'])

    assert all(type(vary) is bool for vary in r1_vary + minf_vary + m0_vary), \
        "All vary parameters should be boolean values."

    return const_dict, pars_dict


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

    return {
        'title': data_title,
        'n_points': n_points,
        'time_points': time_points,
        'magnetizations': np.array(data_values)
    }


def print_parameters(const_dict, pars_dict):
    print(f"Number of sites: {const_dict['n_sites']}")
    for i in range(const_dict['n_sites']):
        print(f"Site {i+1}:")
        print(f"  R1 relaxation rate: {pars_dict['r1_guess'][i]:.4f}, "
              f"vary: {pars_dict['r1_vary'][i]}")
        print(f"  M_inf: {pars_dict['minf_guess'][i]}, "
              f"vary: {pars_dict['minf_vary'][i]}")
        print(f"  M_0: {pars_dict['m0_guess'][i]}, "
              f"vary: {pars_dict['m0_vary'][i]}")
        print()

    print(f"Number of processes: {const_dict['n_procs']}")
    for i in range(const_dict['n_procs']):
        print(f"Process {i+1}: k = {pars_dict['k_guess'][i]}, "
              f"vary = {pars_dict['k_vary'][i]}")
        print("Process matrix:")
        print(const_dict['matrices'][i])
        print()


def model_magnetization(params, const_dict, time_points):
    """
    Compute the model magnetization for all time points and sites.

    params: lmfit.Parameters object
    const_dict: mechanism dictionary
    time_points: list of time points to calculate magnetization for

    Returns: 2D array of calculated magnetizations.
    """

    n_sites = const_dict['n_sites']
    n_procs = const_dict['n_procs']
    time_points = np.array(time_points)
    n_points = len(time_points)

    # Unpack parameters from lmfit Parameters object
    r1 = np.array([params[f"r1_{i}"] for i in range(n_sites)])
    minf = np.array([params[f"minf_{i}"] for i in range(n_sites)])
    m0minf = np.array([params[f"m0_{i}"] for i in range(n_sites)])-minf
    rates = np.array([params[f"rate_{i}"] for i in range(n_procs)])

    # Build exchange matrix for each process and sum
    exch_mat = np.zeros((n_sites, n_sites))
    for i, mat in enumerate(const_dict['matrices']):
        exch_mat += mat * rates[i]

    # Add relaxation rates to diagonal
    for i in range(n_sites):
        exch_mat[i, i] += r1[i]
        for j in range(n_sites):
            if j != i:
                exch_mat[i, i] -= exch_mat[i, j]

    # Eigen decomposition
    eigvals, eigvecs = np.linalg.eig(exch_mat)
    eigvecs_inv = np.linalg.inv(eigvecs)

    # Calculate magnetization for each time point and site
    mags = np.zeros((n_points, n_sites))
    for i, t in enumerate(time_points):
        exp_diag = np.diag(np.exp(-eigvals * t))
        mags[i, :] = eigvecs @ exp_diag @ eigvecs_inv @ m0minf + minf

    return mags


def sir_residuals(params, const_dict, data_dict):
    time_points = data_dict['time_points']
    mags = model_magnetization(params, const_dict, time_points).flatten()
    obs = np.array(data_dict['magnetizations']).flatten()
    return obs - mags


def do_fit(const_dict, pars_dict, data_dict):
    """
    Perform the nonlinear least squares fit using lmfit.
    """

    n_sites = const_dict['n_sites']
    n_procs = const_dict['n_procs']

    # Build lmfit Parameters object
    params = Parameters()
    for i in range(n_sites):
        params.add(f"r1_{i}", value=pars_dict['r1_guess'][i],
                   vary=pars_dict['r1_vary'][i],
                   max=100*pars_dict['r1_guess'][i], min=0)
        params.add(f"minf_{i}", value=pars_dict['minf_guess'][i],
                   vary=pars_dict['minf_vary'][i], max=1.2,
                   min=0.88)
        params.add(f"m0_{i}", value=pars_dict['m0_guess'][i],
                   vary=pars_dict['m0_vary'][i],
                   max=pars_dict['minf_guess'][i],
                   min=-pars_dict['minf_guess'][i])
    for i in range(n_procs):
        params.add(f"rate_{i}", value=pars_dict['k_guess'][i],
                   vary=pars_dict['k_vary'][i],
                   max=100*pars_dict['k_guess'][i], min=0)

    # Run minimization
    result = minimize(sir_residuals, params, args=(const_dict, data_dict),
                      method='leastsq')  # Use Levenberg-Marquardt minimization

    # Print fit report
    report_fit(result)

    # Update pars_dict with fitted values
    for i in range(n_sites):
        pars_dict['r1_guess'][i] = result.params[f"r1_{i}"].value
        pars_dict['minf_guess'][i] = result.params[f"minf_{i}"].value
        pars_dict['m0_guess'][i] = result.params[f"m0_{i}"].value
    for i in range(n_procs):
        pars_dict['k_guess'][i] = result.params[f"rate_{i}"].value

    data_dict['calculated_magnetizations'] = model_magnetization(
        result.params, const_dict, data_dict['time_points'])

    return result


def calc_mags(const_dict, pars_dict, time_points):
    """
    Calculate the magnetizations for all time points and sites.
    Returns a 2D array of calculated magnetizations.
    """

    params = Parameters()
    for i in range(const_dict['n_sites']):
        params.add(f"r1_{i}", value=pars_dict['r1_guess'][i],
                   vary=pars_dict['r1_vary'][i])
        params.add(f"minf_{i}", value=pars_dict['minf_guess'][i],
                   vary=pars_dict['minf_vary'][i])
        params.add(f"m0_{i}", value=pars_dict['m0_guess'][i],
                   vary=pars_dict['m0_vary'][i])
    for i in range(const_dict['n_procs']):
        params.add(f"rate_{i}", value=pars_dict['k_guess'][i],
                   vary=pars_dict['k_vary'][i])

    return model_magnetization(params, const_dict, time_points)


def write_to_csv(filename, const_dict, pars_dict, data_dict):
    """
    Write calculated, observed, and difference values, plus a smooth curve,
    to a CSV file.
    """

    with open(filename, "w", newline="") as csvfile:
        # writer = csv.writer(csvfile, delimiter='')
        # writer.writerow(["Calculated, Observed, and Difference Values"])
        csvfile.write("Calculated, Observed, and Difference Values\n")

        calc = data_dict['calculated_magnetizations']
        for i, time in enumerate(data_dict['time_points']):
            row = f"{time:8.5f}, "

            # Calculated
            for j in range(const_dict['n_sites']):
                row += f" {calc[i, j]:8.4f}, "
            # Observed
            for j in range(const_dict['n_sites']):
                row += f" {data_dict['magnetizations'][i, j]:8.4f}, "
            # Differences
            for j in range(const_dict['n_sites']):
                diff = calc[i, j] - data_dict['magnetizations'][i, j]
                row += f" {diff:8.4f}, "
            csvfile.write(row + '\n')
            # writer.writerow(row)

        # Smooth curve
        csvfile.write("\nCalculated Smooth Curve\n")

        time_0 = 0.0
        time_f = data_dict['time_points'][-1]
        times_smooth = np.linspace(time_0, time_f, 101)[np.newaxis]
        mags_smooth = calc_mags(const_dict, pars_dict, times_smooth.flatten())

        smooth_data = np.concatenate((times_smooth.T, mags_smooth), axis=1)
        for smooth_row in smooth_data:
            csvfile.write(f"{smooth_row[0]:8.5f}, ")
            for val in smooth_row[1:]:
                csvfile.write(f" {val:8.4f}, ")
            csvfile.write('\n')


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
    print(f"sirfit {VERSION}")

    data_filename = f"{args.filename}.dat"
    mech_filename = f"{args.filename}.mch"
    out_filename = f"{args.filename}.out"
    csv_filename = f"{args.filename}.csv"

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
        with open(out_filename, 'w') as outfile:
            outfile.write(f"sirfit {VERSION}")
            outfile.write(f"Mechanism file: {mech_filename}")
            outfile.write(f"Data file: {data_filename}")

            # print(f"Output file: {out_filename}")

    except Exception as e:
        print(f"Error creating output file: {e}")
        return 1

    # Read mechanism and data files
    const_dict, pars_dict = read_mch_file(mech_filename,
                                          specify_vary=args.specify_vary)
    data_dict = read_data_file(data_filename)
    print_parameters(const_dict, pars_dict)

    # Run fitting
    do_fit(const_dict, pars_dict, data_dict)

    # Output results
    print("Fitting complete. Results:")
    print_parameters(const_dict, pars_dict)
    write_to_csv(csv_filename, const_dict, pars_dict, data_dict)


if __name__ == "__main__":
    main()
