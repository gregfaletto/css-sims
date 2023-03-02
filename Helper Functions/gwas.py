import h5py
import numpy as np
# h5py.run_tests()

filename = "4.hdf5"

num_snps = 50000

with h5py.File(filename, "r") as f:
    # List all groups
    print("Keys: %s" % f.keys())
    accessions_key = list(f.keys())[0]
    positions_key = list(f.keys())[1]
    snps_key = list(f.keys())[2]

    # Get the data
    accessions_ref = f[accessions_key]
    print("accessions_ref:")
    print(accessions_ref)
    data_accessions = list(accessions_ref)
    print(data_accessions)
    print(dir(data_accessions))
    print(data_accessions.__class__)
    print(len(data_accessions))
    print(format(data_accessions))



    print("Reading positions_ref...")
    positions_ref = f[positions_key]
    print("positions_ref:")
    print(positions_ref)
    print("Reading positions[0:10]...")
    data_positions = list(positions_ref[0:10])
    print("Done!")
    print(dir(data_positions))
    print(data_positions)
    print(data_positions.__class__)
    print(len(data_positions))
    print(format(data_positions))

    print("Reading snps_ref...")
    snps_ref = f[snps_key]
    print("snps_ref:")
    print(snps_ref)
    print("Reading snps_ref[0:num_snps, :]...")
    data_snps = list(snps_ref[0:num_snps, :])
    # data_snps = list(snps_ref[(10709466 - num_snps):(10709466 + 1), :])
    # data_snps = list(snps_ref[(500000):(500000 + num_snps), :])
    print("Done!")
    print("dir(data_snps)")
    print(dir(data_snps))
    # print(data_snps)
    print("data_snps.__class__")
    print(data_snps.__class__)
    print("len(data_snps)")
    print(len(data_snps))
    # print("format(data_snps)")
    # print(format(data_snps))
    # print(data_snps[0])

    # Save snps to disk as .csv
    np.savetxt("snps50000.csv", data_snps, delimiter=",")

    # After you are done
    f.close()


