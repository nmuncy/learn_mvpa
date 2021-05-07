# %%
"""
Notes

1) Script will generate named sbatch subprocesses named for each sub-step.
    Output for each sbatch step is written to work_dir and prepeneded
    with "sbatch_writeOut"

2) Pre-processing steps include copying data to derivatives, distortion correcting,
    volume registration, normalization, generating intersection masks, blurring, and
    scaling the data.

3) Written in Python 3.8, has afni and c3d dependencies.
"""

import json
import os
import subprocess
import fnmatch
import math
from argparse import ArgumentParser
from gp_step0_dcm2nii import func_sbatch


# make a list of all EPI scans of phase X
def func_epi_list(phase, h_dir):
    h_list = [
        x.split("+")[0]
        for x in os.listdir(h_dir)
        if fnmatch.fnmatch(x, f"run-*_{phase}+orig.HEAD")
    ]
    h_list.sort()
    return h_list


def flatten_list(list_2d):
    flat_list = []
    for element in list_2d:
        if type(element) is list:
            for item in element:
                flat_list.append(item)
        else:
            flat_list.append(element)
    return flat_list


# receive arguments
def func_argparser():
    parser = ArgumentParser("Receive Bash args from wrapper")
    parser.add_argument("h_sub", help="Subject ID")
    parser.add_argument("h_ses", help="Session")
    parser.add_argument("h_par", help="Parent Directory")
    parser.add_argument("h_phl", nargs="+", help="Phase List")
    return parser


# %%
# def func_preproc(data_dir, work_dir, subj, sess, phase_list):

"""
Step 1: Copy data into work_dir

1) Get func, anat, fmap data. Rename appropriately.

2) To account for different num of fmaps, will
    produce AP, PA fmap per run.
"""

# For testing
subj = "sub-005"
sess = "ses-S1"
phase_list = ["loc", "Study"]

par_dir = "/scratch/madlab/nate_vCAT"
data_dir = os.path.join(par_dir, "dset", subj, sess)
work_dir = os.path.join(par_dir, "derivatives", subj, sess)

if not os.path.exists(work_dir):
    os.makedirs(work_dir)

# %%
# Start
subj_num = subj.split("-")[1]

# struct
if not os.path.exists(os.path.join(work_dir, "struct+orig.HEAD")):
    h_cmd = f"""
        module load afni-20.2.06
        3dcopy \
            {data_dir}/anat/{subj}_{sess}_T1w.nii.gz \
            {work_dir}/struct+orig
    """
    h_copy = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    h_copy.wait()

# epi - keep runs sorted by phase
for phase in phase_list:

    epiNii_list = [
        epi
        for epi in os.listdir(os.path.join(data_dir, "func"))
        if fnmatch.fnmatch(epi, f"*task-{phase}*.nii.gz")
    ]
    epiNii_list.sort()

    for h_count, h_file in enumerate(epiNii_list):
        run_count = 1 + h_count
        if not os.path.exists(f"{work_dir}/run-{run_count}_{phase}+orig.HEAD"):
            h_cmd = f"""
                module load afni-20.2.06
                3dcopy \
                    {data_dir}/func/{h_file} \
                    {work_dir}/run-{run_count}_{phase}+orig
            """
            h_copy = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
            h_copy.wait()

# %%
# fmap
fmap_list = [
    x
    for x in os.listdir(os.path.join(data_dir, "fmap"))
    if fnmatch.fnmatch(x, "*.nii.gz")
]

for fmap in fmap_list:
    h_dir = fmap.split("-")[-1].split("_")[0]
    if not os.path.exists(f"{work_dir}/blip_{h_dir}+orig.HEAD"):
        h_cmd = f"""
            module load afni-20.2.06
            3dcopy \
                {data_dir}/fmap/{fmap} \
                {work_dir}/blip_{h_dir}
        """
        h_copy = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        h_copy.wait()


# %%
"""
Step 2: Detect outliers voxels, blip correct

1) Determine the proportion of voxels per volume that have outlier signal.
    Censor volumes that exceed limit.

2) Correct for signal fallout using fmap. This approach is taken from
    afni_proc. It uses the fmap to "unwarp" the run epi.
"""

for phase in phase_list:

    epi_list = func_epi_list(phase, work_dir)
    for run in epi_list:

        # calc tr length
        h_cmd = f"module load afni-20.2.06 \n cd {work_dir} \n 3dinfo -tr {run}+orig"
        h_tr = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        h_len_tr = h_tr.communicate()[0]
        len_tr = h_len_tr.decode("utf-8").strip()

        # calc number of volumes
        h_cmd = (
            f"module load afni-20.2.06 \n cd {work_dir} \n 3dinfo -ntimes {run}+orig"
        )
        h_nvol = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        h_nvol_out = h_nvol.communicate()[0]
        num_tr = h_nvol_out.decode("utf-8").strip()

        # determine polort argument
        pol = 1 + math.ceil((float(len_tr) * float(num_tr)) / 150)

        # determine percentage outliers for e/volume
        if not os.path.exists(os.path.join(work_dir, f"out.cen.{run}.1D")):
            h_cmd = f"""
                cd {work_dir}
                3dToutcount -automask -fraction -polort {pol} \
                    -legendre {run}+orig > outcount.{run}.1D
                1deval -a outcount.{run}.1D -expr '1-step(a-0.05)' > out.cen.{run}.1D
            """
            func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}out", work_dir)

# %%
# Get all epi to be corrected w/same fmap
epiAll_list = [
    x.split("+")[0]
    for x in os.listdir(os.path.join(work_dir))
    if fnmatch.fnmatch(x, "run-*HEAD")
]

# new fmap correct
#   -f = same direction as epi run, will
#       become "Forward".
#   Runs several 10s of minutes, crashes
#       when subitted through func_sbatch.
if not os.path.exists(os.path.join(work_dir, "blip_WARP+orig.HEAD")):
    h_cmd = f"""
        module load afni-20.2.06
        cd {work_dir}

        unWarpEPI.py \
            -f blip_PA+orig \
            -r blip_AP+orig \
            -d '{",".join(epiAll_list)}' \
            -a struct+orig \
            -s fmap

        3dcopy \
            unWarpOutput_fmap/03_fmap_MidWarped_Forward_WARP.nii.gz \
            blip_WARP
    """
    print(h_cmd)
    h_fmap = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    h_fmap.wait()

# %%
# copy fmap corrected files to work_dir
#   from unWarpOutput_foo
# 06* is corrected file
fmapCorr_list = [
    x
    for x in os.listdir(os.path.join(work_dir, "unWarpOutput_fmap"))
    if fnmatch.fnmatch(x, "06*.nii.gz")
]

for fmap_file in fmapCorr_list:
    h_run = fmap_file.split("_")[2]
    h_phase = fmap_file.split("_")[3]
    out_file = os.path.join(work_dir, f"{h_run}_{h_phase}_blip")
    if not os.path.exists(f"{out_file}+orig.HEAD"):
        h_cmd = f"""
            module load afni-20.2.06
            3dcopy \
                {work_dir}/unWarpOutput_fmap/{fmap_file} \
                {out_file}
        """
        h_copy = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        h_copy.wait()


# %%
"""
Step 3: Make volreg base

1) The volreg base (epi_vrBase) is the single volume in entire experiment
    with the smallest number of outlier volumes.

Note: If data is not time shifted, do so before determining volreg_base.
    Also, it was easier to write a small bash script due to the multiple
    afni commands.
"""

# combine all outcount.* into one master file (outcount_all.1D)
#      1D to rule them all, and in the darkness bind them
#
# make sure things are in order
out_dict = {}
for phase in phase_list:
    h_list = [
        os.path.join(work_dir, x)
        for x in os.listdir(work_dir)
        if fnmatch.fnmatch(x, f"outcount.*{phase}.1D")
    ]
    h_list.sort()
    out_dict[phase] = h_list
out_list = flatten_list(list(out_dict.values()))

out_all = os.path.join(work_dir, "outcount_all.1D")
with open(out_all, "w") as outfile:
    for out_file in out_list:
        with open(out_file) as infile:
            outfile.write(infile.read())

# %%
# all epi runs in experiment
scan_dict = {}
for phase in phase_list:
    h_list = [
        x.split(".")[0]
        for x in os.listdir(work_dir)
        if fnmatch.fnmatch(x, f"run-*{phase}_blip+orig.HEAD")
    ]
    h_list.sort()
    scan_dict[phase] = h_list
scan_list = flatten_list(list(scan_dict.values()))

# list of volume nums
num_vols = []
for scan in scan_list:
    h_cmd = f"module load afni-20.2.06 \n cd {work_dir} \n 3dinfo -ntimes {scan}"
    h_nvol = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    h_nvol_out = h_nvol.communicate()[0]
    num_tr = h_nvol_out.decode("utf-8").strip()
    num_vols.append(num_tr)

# determine index of min
h_cmd = f"""
    module load afni-20.2.06
    3dTstat -argmin -prefix - {work_dir}/outcount_all.1D\\'
"""
h_mind = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
h_ind = h_mind.communicate()[0]
ind_min = h_ind.decode("utf-8").strip()

# determine min volume
h_cmd = f"""
    module load afni-20.2.06
    1d_tool.py \
        -set_run_lengths {" ".join(num_vols)} \
        -index_to_run_tr {ind_min}
"""
h_minV = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
h_min = h_minV.communicate()[0]
min_runVol = h_min.decode("utf-8").strip().split()

# determine run, volume
min_run = scan_list[int(min_runVol[0])]
min_vol = int(min_runVol[1])

# make epi volreg base
h_cmd = f"""
    module load afni-20.2.06
    cd {work_dir}
    3dbucket -prefix epi_vrBase {min_run}"[{min_vol}]"
"""
h_vrb = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
h_vrb.wait()


# %%
"""
Step 4: Calc, Perfrom normalization

1) This step will perform the rigid alignments of T1-EPI (A)
    EPI-EPI base volume (B), and non-linear diffeomorphic of
    T1-Template (C). Note blip distortion map (D).

2) It will then concatenate these warp matrices, and warp EPI data from
    raw/native space to template space via W=A'+B+C+D. Thus, only one
    interpolation of the epi data occurs.

3) Will also censor volumes that have outlier, to not bias scaling
"""

# Calculate T1-EPI rigid, T1-Template diffeo
atlas_dir = "/home/data/madlab/atlases/vold2_mni"
h_cmd = f"""
    cd {work_dir}

    align_epi_anat.py \
        -anat2epi \
        -anat struct+orig \
        -save_skullstrip \
        -suffix _al_junk \
        -epi epi_vrBase+orig \
        -epi_base 0 \
        -epi_strip 3dAutomask \
        -cost lpc+ZZ \
        -giant_move \
        -check_flip \
        -volreg off \
        -tshift off

    auto_warp.py -base {os.path.join(atlas_dir, "vold2_mni_brain+tlrc")} \
        -input struct_ns+orig -skull_strip_input no

    3dbucket -DAFNI_NIFTI_VIEW=tlrc -prefix struct_ns awpy/struct_ns.aw.nii*
    cp awpy/anat.un.aff.Xat.1D .
    cp awpy/anat.un.aff.qw_WARP.nii .
"""
if not os.path.exists(os.path.join(work_dir, "struct_ns+tlrc.HEAD")):
    func_sbatch(h_cmd, 1, 4, 4, f"{subj_num}dif", work_dir)

for phase in phase_list:
    epi_list = func_epi_list(phase, work_dir)
    for run in epi_list:

        # determine input
        if blip_tog == 1:
            h_in = f"{run}_blip+orig"
        else:
            h_in = f"{run}+orig"

        # Calculate volreg for e/run
        h_cmd = f"""
            cd {work_dir}
            3dvolreg -verbose \
                -zpad 1 \
                -base epi_vrBase+orig \
                -1Dfile dfile.{run}.1D \
                -prefix {run}_volreg \
                -cubic \
                -1Dmatrix_save mat.{run}.vr.aff12.1D \
                {h_in}
        """
        if not os.path.exists(os.path.join(work_dir, f"mat.{run}.vr.aff12.1D")):
            func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}vre", work_dir)

        # calc voxel dimension i
        h_cmd = f"module load afni-20.2.06 \n cd {work_dir} \n 3dinfo -di {run}+orig"
        h_gs = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        h_gs_out = h_gs.communicate()[0]
        grid_size = h_gs_out.decode("utf-8").strip()

        # determine input
        if blip_tog == 1:
            h_nwarp = f"anat.un.aff.qw_WARP.nii mat.{run}.warp.aff12.1D blip_{run.split('_')[0]}_For_WARP+orig"
        else:
            h_nwarp = f"anat.un.aff.qw_WARP.nii mat.{run}.warp.aff12.1D"

        # concat warp calcs, apply to native EPI.
        #   Make a mask
        h_cmd = f"""
            cd {work_dir}

            cat_matvec -ONELINE \
                anat.un.aff.Xat.1D \
                struct_al_junk_mat.aff12.1D -I \
                mat.{run}.vr.aff12.1D > mat.{run}.warp.aff12.1D

            3dNwarpApply -master struct_ns+tlrc \
                -dxyz {grid_size} \
                -source {run}+orig \
                -nwarp '{h_nwarp}' \
                -prefix {run}_warp

            3dcalc -overwrite -a {run}_blip+orig -expr 1 -prefix tmp_{run}_mask

            3dNwarpApply -master struct_ns+tlrc \
                -dxyz {grid_size} \
                -source tmp_{run}_mask+orig \
                -nwarp 'anat.un.aff.qw_WARP.nii mat.{run}.warp.aff12.1D' \
                -interp cubic \
                -ainterp NN -quiet \
                -prefix {run}_mask_warped

            3dTstat -min -prefix tmp_{run}_min {run}_mask_warped+tlrc
        """
        if not os.path.exists(os.path.join(work_dir, f"{run}_warp+tlrc.HEAD")):
            func_sbatch(h_cmd, 2, 4, 4, f"{subj_num}war", work_dir)

# Determine minimum value, make mask
#   beware the jabberwocky i.e. expanding braces in 3dMean
for phase in phase_list:
    epi_list = func_epi_list(phase, work_dir)
    h_cmd = f"""
        cd {work_dir}
        3dMean -datum short -prefix tmp_mean tmp_run-{{1..{len(epi_list)}}}_{phase}_min+tlrc
        3dcalc -a tmp_mean+tlrc -expr 'step(a-0.999)' -prefix {phase}_minVal_mask
    """
    if not os.path.exists(os.path.join(work_dir, f"{phase}_minVal_mask+tlrc.HEAD")):
        func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}min", work_dir)

    # make clean data
    epi_mask = os.path.join(work_dir, f"{phase}_minVal_mask+tlrc")
    for run in epi_list:
        h_warp = os.path.join(work_dir, f"{run}_warp+tlrc")
        h_clean = os.path.join(work_dir, f"{run}_volreg_clean")
        h_cmd = f"3dcalc -a {h_warp} -b {epi_mask} -expr 'a*b' -prefix {h_clean}"
        if not os.path.exists(h_clean + "+tlrc.HEAD"):
            func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}cle", work_dir)

# %%
"""
Step 5: Blur, Make masks

1) Blur epi data - 1.5 * voxel dim, rounded up to nearest int. FWHM.

2) Make a union mask, where sufficient signal exists for both T1w
    and T2*w at each voxel for analyses. Incorporated at the
    group-level analysis.

3) Make tissue masks. The WM mask will be used later to derive
    nuissance regressors for the REML.
    Note: this references some custom masks, and is based in
        atropos rather than in AFNIs tiss seg protocol.
"""

# Blur
for phase in phase_list:
    epi_list = func_epi_list(phase, work_dir)
    for run in epi_list:

        # calc voxel dim i
        h_cmd = f"module load afni-20.2.06 \n cd {work_dir} \n 3dinfo -di {run}+orig"
        h_gs = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        h_gs_out = h_gs.communicate()[0]
        grid_size = h_gs_out.decode("utf-8").strip()
        blur_size = math.ceil(1.5 * float(grid_size))

        # do blur
        h_bin = os.path.join(work_dir, f"{run}_volreg_clean+tlrc")
        h_bout = os.path.join(work_dir, f"{run}_blur")
        h_cmd = f"3dmerge -1blur_fwhm {blur_size} -doall -prefix {h_bout} {h_bin}"
        if not os.path.exists(f"{h_bout}+tlrc.HEAD"):
            func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}blur", work_dir)

# Make EPI-T1 union mask (mask_epi_anat)
run_list = [
    x.split(".")[0]
    for x in os.listdir(work_dir)
    if fnmatch.fnmatch(x, "*blur+tlrc.HEAD")
]

if not os.path.exists(os.path.join(work_dir, "mask_epi_anat+tlrc.HEAD")):

    for i in run_list:
        h_mout = os.path.join(work_dir, f"tmp_mask.{i}")
        h_min = os.path.join(work_dir, i)
        if not os.path.exists(h_mout + ".HEAD"):
            h_cmd = f"3dAutomask -prefix {h_mout} {h_min}"
            func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}mau", work_dir)

    h_cmd = f"""
        cd {work_dir}
        3dmask_tool -inputs tmp_mask.*+tlrc.HEAD -union -prefix tmp_mask_allRuns
        3dresample -master tmp_mask_allRuns+tlrc -input struct_ns+tlrc \
            -prefix tmp_anat_resamp
        3dmask_tool -dilate_input 5 -5 -fill_holes -input tmp_anat_resamp+tlrc \
            -prefix tmp_mask_struct
        3dmask_tool -input tmp_mask_allRuns+tlrc tmp_mask_struct+tlrc -inter \
            -prefix mask_epi_anat
        3dABoverlap -no_automask tmp_mask_allRuns+tlrc tmp_mask_struct+tlrc | \
            tee out.mask_ae_overlap.txt
    """
    func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}uni", work_dir)

# Make tissue-class masks
#   I like Atropos better than AFNI's way, so use those priors
atropos_dict = {1: "CSF", 2: "GMc", 3: "WM", 4: "GMs"}
atropos_dir = os.path.join(atlas_dir, "priors_ACT")
h_tcin = os.path.join(work_dir, run_list[0])

for i in atropos_dict:
    h_tiss = atropos_dict[i]
    if not os.path.exists(
        os.path.join(work_dir, f"final_mask_{h_tiss}_eroded+tlrc.HEAD")
    ):
        h_cmd = f"""
            module load c3d-1.0.0-gcc-8.2.0
            cd {work_dir}

            c3d {atropos_dir}/Prior{i}.nii.gz -thresh 0.3 1 1 0 \
                -o tmp_{h_tiss}_bin.nii.gz
            3dresample -master {h_tcin} -rmode NN -input tmp_{h_tiss}_bin.nii.gz \
                -prefix final_mask_{h_tiss}+tlrc
            3dmask_tool -input tmp_{h_tiss}_bin.nii.gz -dilate_input -1 \
                -prefix tmp_mask_{h_tiss}_eroded
            3dresample -master {h_tcin} -rmode NN -input tmp_mask_{h_tiss}_eroded+orig \
                -prefix final_mask_{h_tiss}_eroded
        """
        func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}atr", work_dir)

# %%
"""
Step 6: Scale data

1) Data is scaled by mean signal
"""

for phase in phase_list:
    epi_list = func_epi_list(phase, work_dir)
    for run in epi_list:
        if not os.path.exists(os.path.join(work_dir, f"{run}_scale+tlrc.HEAD")):
            h_cmd = f"""
                cd {work_dir}
                3dTstat -prefix tmp_tstat_{run} {run}_blur+tlrc
                3dcalc -a {run}_blur+tlrc \
                    -b tmp_tstat_{run}+tlrc \
                    -c {phase}_minVal_mask+tlrc \
                    -expr 'c * min(200, a/b*100)*step(a)*step(b)' \
                    -prefix {run}_scale
            """
            func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}scale", work_dir)


# %%
# def main():

#     args = func_argparser().parse_args()
#     h_data_dir = os.path.join(args.h_par, "dset", args.h_sub, args.h_ses)
#     h_work_dir = os.path.join(args.h_par, "derivatives", args.h_sub, args.h_ses)

#     if not os.path.exists(h_work_dir):
#         os.makedirs(h_work_dir)

#     # print(h_data_dir, h_work_dir, args.h_sub, args.h_ses, args.h_phl)
#     func_preproc(h_data_dir, h_work_dir, args.h_sub, args.h_ses, args.h_phl)


# if __name__ == "__main__":
#     main()

# %%
