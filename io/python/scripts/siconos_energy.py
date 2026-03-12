#!/usr/bin/env @Python_EXECUTABLE@

import sys
import numpy as np
import h5py
import os

import getopt
import siconos.plot_config as sicoplot

plt, enable_plot = sicoplot.choose_backend()


def get_data_energy(f, idx):
    x = f["data"]["energy_work"][:, idx]
    return x


def plot_min_max_cf_work(filename):
    data = h5py.File(filename, "r")

    def get_data_cf_work(f, idx):
        x = f["data"]["cf_work"][:, idx]
        return x

    # x = data["data"]["cf_work"][:]
    # print('x', x)
    times_e = get_data_energy(data, 0)

    times = get_data_cf_work(data, 0)
    if times.shape[0] == 0:
        print("no data on in /data/energy_work")
        exit(0)
    # print('times', times)
    times_set = set(times.tolist())
    # print('times_set', times_set)
    times_list = np.sort(list(times_set))
    # print('times_list', times_list)

    cf_work_normal = get_data_cf_work(data, 2)
    cf_work_tangent = get_data_cf_work(data, 3)
    # print('cf_work_normal',  cf_work_normal)

    id_t_old = 0
    # list_id = []

    cf_work_min_max = []
    for t in times_e:
        cf_work_min_max.append([0.0, 0.0, 0.0, 0.0])

    sum = False
    print("start computation of min and max of cf_work. may take a bit of time...")
    for t in times_list:
        # print('t', t)

        id_t = max(0, np.searchsorted(times, t, side="right"))
        id_t_abs = max(0, np.searchsorted(times_e, t, side="right") - 1)

        # print(id_t, id_t_abs)
        # if id_t_abs > len(cf_work_min_max):
        #     continue
        # if id_t_abs > len(times_e):
        #     continue

        # print(len(cf_work_min_max))
        # print(times_e[id_t_abs])

        id_start = id_t_old
        id_end = id_t
        # print('id_t range ', id_start, id_end)

        cf_work_n = cf_work_normal[id_start:id_end]
        cf_work_t = cf_work_tangent[id_start:id_end]

        # print( "cf_work_n is: ", cf_work_n, cf_work_n.shape )
        # print( "cf_work_t is: ", cf_work_t, cf_work_n.shape )
        if cf_work_n.shape[0] > 0:
            if sum is True:
                cf_work_min_max[id_t_abs] = [
                    np.sum(cf_work_n),
                    np.sum(cf_work_n),
                    np.sum(cf_work_t),
                    np.sum(cf_work_t),
                ]
            else:
                cf_work_min_max[id_t_abs] = [
                    np.max(np.where(cf_work_n > 0, cf_work_n, 0.0)),
                    np.min(np.where(cf_work_n < 0, cf_work_n, 0.0)),
                    np.max(np.where(cf_work_t > 0, cf_work_t, 0.0)),
                    np.min(np.where(cf_work_t < 0, cf_work_t, 0.0)),
                ]
        # print('cf_work_min_max', cf_work_min_max[id_t_abs])
        id_t_old = id_t

    cf_work_min_max = np.array(cf_work_min_max)

    basename = os.path.basename(filename)
    nsl = ""
    if "Fremond" in basename:
        nsl = "(Fremond)"
    if "Newton" in basename:
        nsl = "(Newton)"

    plt.figure(num=4, figsize=(10, 10))
    if sum is True:
        plt.subplot(211)
        plt.grid(visible=True)
        plt.plot(cf_work_min_max[:, 0], label="normal work sum" + nsl)
        plt.legend()
        plt.subplot(212)
        plt.plot(cf_work_min_max[:, 2], label="tangent work sum" + nsl)
        plt.legend()
        plt.xlabel("time")
        plt.ylabel("Energy/Work")
        plt.grid(visible=True)
    else:
        plt.suptitle("Siconos min and max contact work plots from hdf5 file")
        plt.subplot(411)
        plt.grid(visible=True)
        plt.plot(times_e, cf_work_min_max[:, 0], label="normal work max" + nsl)
        plt.legend()
        plt.subplot(412)
        plt.plot(times_e, cf_work_min_max[:, 1], label="normal work min" + nsl)
        plt.legend()
        plt.grid(visible=True)
        plt.subplot(413)
        plt.grid(visible=True)
        plt.plot(times_e, cf_work_min_max[:, 2], label="tangent work max" + nsl)
        plt.legend()
        plt.subplot(414)
        plt.plot(times_e, cf_work_min_max[:, 3], label="tangent work min" + nsl)
        plt.legend()
        plt.grid(visible=True)
        plt.xlabel("time")
        plt.ylabel("Energy/Work")
        plt.savefig("./figures_energy/" + basename + "_work_min_max.png")

    print("max tangent work max", np.max(cf_work_min_max[:, 2]))
    print("min tangent work min", np.min(cf_work_min_max[:, 2]))


def plot_energy(filename):
    data = h5py.File(filename, "r")

    time = get_data_energy(data, 0)
    if time.shape[0] == 0:
        print("no data on in /data/energy_work")
        exit(0)
    kinetic = get_data_energy(data, 1)

    force_work = np.cumsum(get_data_energy(data, 2))
    force_work = force_work  # - force_work[-1]
    normal_contact_work = np.cumsum(get_data_energy(data, 3))
    tangent_contact_work = np.cumsum(get_data_energy(data, 4))

    normal_contact_work_negative = np.cumsum(get_data_energy(data, 7))
    tangent_contact_work_negative = np.cumsum(get_data_energy(data, 8))

    basename = os.path.basename(filename)
    nsl = ""
    if "Fremond" in basename:
        nsl = "(Fremond)"
    if "Newton" in basename:
        nsl = "(Newton)"

    start_index = 1

    plt.figure(num=1, figsize=(20, 10))
    plt.suptitle("Siconos energy/work plots from hdf5 file")
    plt.subplot(121)
    plt.plot(
        time[start_index:],
        force_work[start_index:],
        label="force work (potential energy)" + nsl,
    )

    plt.plot(time[start_index:], kinetic[start_index:], label="kinetic" + nsl)
    plt.plot(
        time[start_index:],
        normal_contact_work[start_index:],
        label="normal contact work" + nsl,
    )
    plt.plot(
        time[start_index:],
        tangent_contact_work[start_index:],
        label="tangent contact work" + nsl,
    )

    plt.plot(
        time[start_index:],
        -force_work[start_index:]
        + kinetic[start_index:]
        - normal_contact_work[start_index:]
        - tangent_contact_work[start_index:],
        label="energy_balance" + nsl,
        linestyle="dotted",
    )
    plt.xlabel("time")
    plt.ylabel("Energy/Work")
    plt.legend()

    plt.grid(visible=True)
    plt.subplot(122)
    plt.plot(
        time[start_index:],
        force_work[start_index:],
        label="force work (potential energy)" + nsl,
    )
    plt.fill_between(
        time[start_index:],
        force_work[start_index:],
        force_work[start_index:] - kinetic[start_index:],
        alpha=0.3,
        label="kinetic energy",
    )
    plt.plot(time[start_index:], force_work[start_index:] - kinetic[start_index:])
    plt.fill_between(
        time[start_index:],
        force_work[start_index:] - kinetic[start_index:],
        force_work[start_index:]
        - kinetic[start_index:]
        + normal_contact_work[start_index:],
        alpha=0.3,
        label="normal contact work",
    )
    plt.plot(
        time[start_index:],
        force_work[start_index:]
        - kinetic[start_index:]
        + normal_contact_work[start_index:],
    )
    plt.fill_between(
        time[start_index:],
        force_work[start_index:]
        - kinetic[start_index:]
        + normal_contact_work[start_index:],
        force_work[start_index:]
        - kinetic[start_index:]
        + normal_contact_work[start_index:]
        + tangent_contact_work[start_index:],
        alpha=0.3,
        label="tangent contact work",
    )

    plt.plot(
        time[start_index:],
        force_work[start_index:]
        - kinetic[start_index:]
        + normal_contact_work[start_index:]
        + tangent_contact_work[start_index:],
    )
    plt.fill_between(
        time[start_index:],
        force_work[start_index:]
        - kinetic[start_index:]
        + normal_contact_work[start_index:]
        + tangent_contact_work[start_index:],
        np.zeros_like(time[start_index:]),
        alpha=0.3,
        label="numerical loss",
    )

    plt.xlabel("time")
    plt.ylabel("Energy/Work")
    plt.legend()

    plt.grid(visible=True)

    if with_save:
        if not os.path.exists("./figures_energy"):
            os.makedirs("./figures_energy")
        plt.savefig("./figures_energy/" + basename + "_energy.png")

    if with_details:
        plt.figure(num=2, figsize=(20, 10))
        plt.suptitle("Siconos energy/work plots details from hdf5 file")
        plt.subplot(121)
        plt.plot(
            time[start_index:],
            normal_contact_work[start_index:],
            label="normal contact work" + nsl,
        )
        plt.plot(
            time[start_index:],
            tangent_contact_work[start_index:],
            label="tangent contact work" + nsl,
        )
        plt.plot(
            time[start_index:],
            normal_contact_work_negative[start_index:],
            label="normal contact work (only negative part)" + nsl,
        )
        plt.plot(
            time[start_index:],
            tangent_contact_work_negative[start_index:],
            label="tangent contact work (only negative part)" + nsl,
        )
        plt.legend()
        plt.xlabel("time")
        plt.ylabel("Energy/Work")
        plt.grid(visible=True)

        plt.subplot(122)
        plt.plot(
            time[start_index:],
            -force_work[start_index:]
            + kinetic[start_index:]
            - normal_contact_work[start_index:]
            - tangent_contact_work[start_index:],
            label="energy_balance" + nsl,
            linestyle="dotted",
        )
        plt.legend()
        plt.grid(visible=True)
        plt.xlabel("time")
        plt.ylabel("Energy/Work")

        if with_save:
            if not os.path.exists("./figures_energy"):
                os.makedirs("./figures_energy")
            plt.savefig("./figures_energy/" + basename + "_energy_balance_detail.png")

        # plt.figure(num=4,figsize=(10,10))
        # plt.plot(time[start_index:], -force_work[start_index:], label='force_work'+nsl)
        # plt.legend()
        # plt.grid(visible=True)


def usage(long=False):
    print()
    print("Usage: {0} [OPTION]... <HDF5>".format(os.path.split(sys.argv[0])[1]))

    print("""[--help] [--details] [--save] """)

    print()


try:
    opts, args = getopt.gnu_getopt(
        sys.argv[1:],
        "",
        ["help", "details", "save", "min-max-cf_work"],
    )

except getopt.GetoptError as err:
    sys.stderr.write("{0}\n".format(str(err)))
    usage()
    exit(2)

with_details = False
with_save = False
with_min_max_cf_work = False
for o, a in opts:
    if o == "--help":
        usage(long=True)
        exit(0)
    if o == "--details":
        with_details = True
    if o == "--save":
        with_save = True
    if o == "--min-max-cf_work":
        with_min_max_cf_work = True

if len(args) < 1:
    usage()
    exit(0)
else:
    filename = args[0]


directory = ""
plot_energy(filename)

if with_min_max_cf_work:
    plot_min_max_cf_work(filename)

if enable_plot:
    plt.show()
else:
    plt.savefig(filename + ".png")
