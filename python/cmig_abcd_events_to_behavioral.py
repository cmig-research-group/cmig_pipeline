#!/usr/bin/env python
"""
cmig_abcd_events_to_behavioral.py

Reads one or more BIDS _events.tsv files (one per run) and writes a
wide-format behavioral metrics TSV.  The task is read from each file
name (task-<label>), and metadata from constant-valued columns in the TSV
(task_version, etc.).
For MID/SST/NBACK it also writes the E-Prime experiment name and versions
and the min/max of the runs' timing QC delays (func/ runs only).
Does not access the original E-Prime file.

Public function
---------------
events_to_behavioral(fname_tsv_list, fname_json=None,
                     outdir='.', outstem=None, run_slots=None,
                     forceflag=False, verbose=True)
Returns (behav_dict, errcode, errmsg).

Created:  05/18/2026 by Don Hagler
Prev mod: 06/25/2026 by Don Hagler
Prev mod: 06/29/2026 by Don Hagler
Prev mod: 07/07/2026 by Don Hagler
Prev mod: 07/08/2026 by Don Hagler
Prev mod: 07/09/2026 by Don Hagler
Prev mod: 07/10/2026 by Don Hagler
Prev mod: 07/13/2026 by Don Hagler
Prev mod: 07/27/2026 by Don Hagler
Prev mod: 07/28/2026 by Don Hagler
Prev mod: 10/01/2026 by Don Hagler
Prev mod: 10/20/2026 by Don Hagler
Prev mod: 10/21/2026 by Don Hagler
Prev mod: 11/01/2026 by Don Hagler
Prev mod: 11/12/2026 by Don Hagler
Last mod: 11/15/2026 by Don Hagler
"""
from __future__ import print_function
import math
import os
import re

import numpy as np
import pandas as pd
from scipy.stats import norm as _norm


# ── release-name mapping (matlab_name → ABCD release column name) ──────────
# Keys that are absent from this dict (SubjID, VisitID, version) are
# dropped from the output TSV.
_RELEASE_NAMES = {
    'MID': {
        'switch_flag': 'mr_y_tfmri__mid__beh__switch_indicator',
        'perform_flag': 'mr_y_tfmri__mid__beh__qc_indicator',
        'nruns': 'mr_y_tfmri__mid__beh__run_count',
        'eprime_behav_only': 'mr_y_tfmri__mid__beh__eprime__behav_indicator',
        'eprime_expt_name': 'mr_y_tfmri__mid__beh__eprime__expt_name',
        'eprime_expt_version': 'mr_y_tfmri__mid__beh__eprime__expt_ver',
        'eprime_runtime_version': 'mr_y_tfmri__mid__beh__eprime__runtime_ver',
        'start_delay_min': 'mr_y_tfmri__mid__beh__start__delay_min',
        'start_delay_max': 'mr_y_tfmri__mid__beh__start__delay_max',
        'start_delay_diff_min': 'mr_y_tfmri__mid__beh__start__delay__diff_min',
        'start_delay_diff_max': 'mr_y_tfmri__mid__beh__start__delay__diff_max',
        'init_delay_diff_min': 'mr_y_tfmri__mid__beh__init__delay__diff_min',
        'init_delay_diff_max': 'mr_y_tfmri__mid__beh__init__delay__diff_max',
        'total_numtrials': 'mr_y_tfmri__mid__beh__trial_count',
        'total_numtrials_run_1': 'mr_y_tfmri__mid__beh__r01__trial_count',
        'total_numtrials_run_2': 'mr_y_tfmri__mid__beh__r02__trial_count',
        'small_reward_total_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__rwrd__sm_count',
        'small_reward_total_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__rwrd__sm_rt',
        'small_reward_total_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__rwrd__sm_std',
        'small_reward_total_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__sm_count',
        'small_reward_total_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__sm_rt',
        'small_reward_total_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__sm_std',
        'small_reward_total_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__sm_count',
        'small_reward_total_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__sm_rt',
        'small_reward_total_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__sm_std',
        'small_reward_pos_feedback_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__rwrd__sm__crct_count',
        'small_reward_pos_feedback_bothruns_rate': 'mr_y_tfmri__mid__beh__trial__rwrd__sm__crct_prop',
        'small_reward_pos_feedback_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__rwrd__sm__crct_rt',
        'small_reward_pos_feedback_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__rwrd__sm__crct_std',
        'small_reward_pos_feedback_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__sm__crct_count',
        'small_reward_pos_feedback_run_1_rate': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__sm__crct_prop',
        'small_reward_pos_feedback_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__sm__crct_rt',
        'small_reward_pos_feedback_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__sm__crct_std',
        'small_reward_pos_feedback_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__sm__crct_count',
        'small_reward_pos_feedback_run_2_rate': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__sm__crct_prop',
        'small_reward_pos_feedback_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__sm__crct_rt',
        'small_reward_pos_feedback_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__sm__crct_std',
        'small_reward_neg_feedback_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__rwrd__sm__incrct_count',
        'small_reward_neg_feedback_bothruns_rate': 'mr_y_tfmri__mid__beh__trial__rwrd__sm__incrct_prop',
        'small_reward_neg_feedback_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__rwrd__sm__incrct_rt',
        'small_reward_neg_feedback_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__rwrd__sm__incrct_std',
        'small_reward_neg_feedback_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__sm__incrct_count',
        'small_reward_neg_feedback_run_1_rate': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__sm__incrct_prop',
        'small_reward_neg_feedback_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__sm__incrct_rt',
        'small_reward_neg_feedback_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__sm__incrct_std',
        'small_reward_neg_feedback_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__sm__incrct_count',
        'small_reward_neg_feedback_run_2_rate': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__sm__incrct_prop',
        'small_reward_neg_feedback_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__sm__incrct_rt',
        'small_reward_neg_feedback_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__sm__incrct_std',
        'large_reward_total_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__rwrd__lg_count',
        'large_reward_total_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__rwrd__lg_rt',
        'large_reward_total_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__rwrd__lg_std',
        'large_reward_total_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__lg_count',
        'large_reward_total_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__lg_rt',
        'large_reward_total_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__lg_std',
        'large_reward_total_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__lg_count',
        'large_reward_total_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__lg_rt',
        'large_reward_total_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__lg_std',
        'large_reward_pos_feedback_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__rwrd__lg__crct_count',
        'large_reward_pos_feedback_bothruns_rate': 'mr_y_tfmri__mid__beh__trial__rwrd__lg__crct_prop',
        'large_reward_pos_feedback_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__rwrd__lg__crct_rt',
        'large_reward_pos_feedback_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__rwrd__lg__crct_std',
        'large_reward_pos_feedback_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__lg__crct_count',
        'large_reward_pos_feedback_run_1_rate': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__lg__crct_prop',
        'large_reward_pos_feedback_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__lg__crct_rt',
        'large_reward_pos_feedback_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__lg__crct_std',
        'large_reward_pos_feedback_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__lg__crct_count',
        'large_reward_pos_feedback_run_2_rate': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__lg__crct_prop',
        'large_reward_pos_feedback_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__lg__crct_rt',
        'large_reward_pos_feedback_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__lg__crct_std',
        'large_reward_neg_feedback_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__rwrd__lg__incrct_count',
        'large_reward_neg_feedback_bothruns_rate': 'mr_y_tfmri__mid__beh__trial__rwrd__lg__incrct_prop',
        'large_reward_neg_feedback_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__rwrd__lg__incrct_rt',
        'large_reward_neg_feedback_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__rwrd__lg__incrct_std',
        'large_reward_neg_feedback_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__lg__incrct_count',
        'large_reward_neg_feedback_run_1_rate': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__lg__incrct_prop',
        'large_reward_neg_feedback_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__lg__incrct_rt',
        'large_reward_neg_feedback_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__lg__incrct_std',
        'large_reward_neg_feedback_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__lg__incrct_count',
        'large_reward_neg_feedback_run_2_rate': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__lg__incrct_prop',
        'large_reward_neg_feedback_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__lg__incrct_rt',
        'large_reward_neg_feedback_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__lg__incrct_std',
        'small_loss_total_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__loss__sm_count',
        'small_loss_total_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__loss__sm_rt',
        'small_loss_total_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__loss__sm_std',
        'small_loss_total_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__loss__sm_count',
        'small_loss_total_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__loss__sm_rt',
        'small_loss_total_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__loss__sm_std',
        'small_loss_total_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__loss__sm_count',
        'small_loss_total_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__loss__sm_rt',
        'small_loss_total_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__loss__sm_std',
        'small_loss_pos_feedback_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__loss__sm__crct_count',
        'small_loss_pos_feedback_bothruns_rate': 'mr_y_tfmri__mid__beh__trial__loss__sm__crct_prop',
        'small_loss_pos_feedback_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__loss__sm__crct_rt',
        'small_loss_pos_feedback_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__loss__sm__crct_std',
        'small_loss_pos_feedback_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__loss__sm__crct_count',
        'small_loss_pos_feedback_run_1_rate': 'mr_y_tfmri__mid__beh__r01__trial__loss__sm__crct_prop',
        'small_loss_pos_feedback_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__loss__sm__crct_rt',
        'small_loss_pos_feedback_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__loss__sm__crct_std',
        'small_loss_pos_feedback_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__loss__sm__crct_count',
        'small_loss_pos_feedback_run_2_rate': 'mr_y_tfmri__mid__beh__r02__trial__loss__sm__crct_prop',
        'small_loss_pos_feedback_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__loss__sm__crct_rt',
        'small_loss_pos_feedback_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__loss__sm__crct_std',
        'small_loss_neg_feedback_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__loss__sm__incrct_count',
        'small_loss_neg_feedback_bothruns_rate': 'mr_y_tfmri__mid__beh__trial__loss__sm__incrct_prop',
        'small_loss_neg_feedback_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__loss__sm__incrct_rt',
        'small_loss_neg_feedback_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__loss__sm__incrct_std',
        'small_loss_neg_feedback_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__loss__sm__incrct_count',
        'small_loss_neg_feedback_run_1_rate': 'mr_y_tfmri__mid__beh__r01__trial__loss__sm__incrct_prop',
        'small_loss_neg_feedback_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__loss__sm__incrct_rt',
        'small_loss_neg_feedback_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__loss__sm__incrct_std',
        'small_loss_neg_feedback_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__loss__sm__incrct_count',
        'small_loss_neg_feedback_run_2_rate': 'mr_y_tfmri__mid__beh__r02__trial__loss__sm__incrct_prop',
        'small_loss_neg_feedback_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__loss__sm__incrct_rt',
        'small_loss_neg_feedback_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__loss__sm__incrct_std',
        'large_loss_total_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__loss__lg_count',
        'large_loss_total_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__loss__lg_rt',
        'large_loss_total_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__loss__lg_std',
        'large_loss_total_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__loss__lg_count',
        'large_loss_total_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__loss__lg_rt',
        'large_loss_total_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__loss__lg_std',
        'large_loss_total_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__loss__lg_count',
        'large_loss_total_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__loss__lg_rt',
        'large_loss_total_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__loss__lg_std',
        'large_loss_pos_feedback_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__loss__lg__crct_count',
        'large_loss_pos_feedback_bothruns_rate': 'mr_y_tfmri__mid__beh__trial__loss__lg__crct_prop',
        'large_loss_pos_feedback_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__loss__lg__crct_rt',
        'large_loss_pos_feedback_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__loss__lg__crct_std',
        'large_loss_pos_feedback_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__loss__lg__crct_count',
        'large_loss_pos_feedback_run_1_rate': 'mr_y_tfmri__mid__beh__r01__trial__loss__lg__crct_prop',
        'large_loss_pos_feedback_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__loss__lg__crct_rt',
        'large_loss_pos_feedback_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__loss__lg__crct_std',
        'large_loss_pos_feedback_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__loss__lg__crct_count',
        'large_loss_pos_feedback_run_2_rate': 'mr_y_tfmri__mid__beh__r02__trial__loss__lg__crct_prop',
        'large_loss_pos_feedback_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__loss__lg__crct_rt',
        'large_loss_pos_feedback_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__loss__lg__crct_std',
        'large_loss_neg_feedback_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__loss__lg__incrct_count',
        'large_loss_neg_feedback_bothruns_rate': 'mr_y_tfmri__mid__beh__trial__loss__lg__incrct_prop',
        'large_loss_neg_feedback_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__loss__lg__incrct_rt',
        'large_loss_neg_feedback_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__loss__lg__incrct_std',
        'large_loss_neg_feedback_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__loss__lg__incrct_count',
        'large_loss_neg_feedback_run_1_rate': 'mr_y_tfmri__mid__beh__r01__trial__loss__lg__incrct_prop',
        'large_loss_neg_feedback_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__loss__lg__incrct_rt',
        'large_loss_neg_feedback_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__loss__lg__incrct_std',
        'large_loss_neg_feedback_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__loss__lg__incrct_count',
        'large_loss_neg_feedback_run_2_rate': 'mr_y_tfmri__mid__beh__r02__trial__loss__lg__incrct_prop',
        'large_loss_neg_feedback_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__loss__lg__incrct_rt',
        'large_loss_neg_feedback_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__loss__lg__incrct_std',
        'neutral_total_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__neut_count',
        'neutral_total_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__neut_rt',
        'neutral_total_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__neut_std',
        'neutral_total_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__neut_count',
        'neutral_total_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__neut_rt',
        'neutral_total_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__neut_std',
        'neutral_total_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__neut_count',
        'neutral_total_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__neut_rt',
        'neutral_total_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__neut_std',
        'neutral_pos_feedback_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__neut__crct_count',
        'neutral_pos_feedback_bothruns_rate': 'mr_y_tfmri__mid__beh__trial__neut__crct_prop',
        'neutral_pos_feedback_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__neut__crct_rt',
        'neutral_pos_feedback_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__neut__crct_std',
        'neutral_pos_feedback_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__neut__crct_count',
        'neutral_pos_feedback_run_1_rate': 'mr_y_tfmri__mid__beh__r01__trial__neut__crct_prop',
        'neutral_pos_feedback_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__neut__crct_rt',
        'neutral_pos_feedback_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__neut__crct_std',
        'neutral_pos_feedback_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__neut__crct_count',
        'neutral_pos_feedback_run_2_rate': 'mr_y_tfmri__mid__beh__r02__trial__neut__crct_prop',
        'neutral_pos_feedback_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__neut__crct_rt',
        'neutral_pos_feedback_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__neut__crct_std',
        'neutral_neg_feedback_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__neut__incrct_count',
        'neutral_neg_feedback_bothruns_rate': 'mr_y_tfmri__mid__beh__trial__neut__incrct_prop',
        'neutral_neg_feedback_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__neut__incrct_rt',
        'neutral_neg_feedback_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__neut__incrct_std',
        'neutral_neg_feedback_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__neut__incrct_count',
        'neutral_neg_feedback_run_1_rate': 'mr_y_tfmri__mid__beh__r01__trial__neut__incrct_prop',
        'neutral_neg_feedback_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__neut__incrct_rt',
        'neutral_neg_feedback_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__neut__incrct_std',
        'neutral_neg_feedback_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__neut__incrct_count',
        'neutral_neg_feedback_run_2_rate': 'mr_y_tfmri__mid__beh__r02__trial__neut__incrct_prop',
        'neutral_neg_feedback_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__neut__incrct_rt',
        'neutral_neg_feedback_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__neut__incrct_std',
        'hiloreward_total_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__rwrd_count',
        'hiloreward_total_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__rwrd_rt',
        'hiloreward_total_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__rwrd_std',
        'hiloreward_total_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__rwrd_count',
        'hiloreward_total_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__rwrd_rt',
        'hiloreward_total_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__rwrd_std',
        'hiloreward_total_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__rwrd_count',
        'hiloreward_total_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__rwrd_rt',
        'hiloreward_total_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__rwrd_std',
        'hiloreward_pos_feedback_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__rwrd__crct_count',
        'hiloreward_pos_feedback_bothruns_rate': 'mr_y_tfmri__mid__beh__trial__rwrd__crct_prop',
        'hiloreward_pos_feedback_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__rwrd__crct_rt',
        'hiloreward_pos_feedback_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__rwrd__crct_std',
        'hiloreward_pos_feedback_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__crct_count',
        'hiloreward_pos_feedback_run_1_rate': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__crct_prop',
        'hiloreward_pos_feedback_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__crct_rt',
        'hiloreward_pos_feedback_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__crct_std',
        'hiloreward_pos_feedback_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__crct_count',
        'hiloreward_pos_feedback_run_2_rate': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__crct_prop',
        'hiloreward_pos_feedback_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__crct_rt',
        'hiloreward_pos_feedback_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__crct_std',
        'hiloreward_neg_feedback_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__rwrd__incrct_count',
        'hiloreward_neg_feedback_bothruns_rate': 'mr_y_tfmri__mid__beh__trial__rwrd__incrct_prop',
        'hiloreward_neg_feedback_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__rwrd__incrct_rt',
        'hiloreward_neg_feedback_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__rwrd__incrct_std',
        'hiloreward_neg_feedback_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__incrct_count',
        'hiloreward_neg_feedback_run_1_rate': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__incrct_prop',
        'hiloreward_neg_feedback_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__incrct_rt',
        'hiloreward_neg_feedback_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__rwrd__incrct_std',
        'hiloreward_neg_feedback_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__incrct_count',
        'hiloreward_neg_feedback_run_2_rate': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__incrct_prop',
        'hiloreward_neg_feedback_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__incrct_rt',
        'hiloreward_neg_feedback_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__rwrd__incrct_std',
        'hiloloss_total_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__loss_count',
        'hiloloss_total_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__loss_rt',
        'hiloloss_total_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__loss_std',
        'hiloloss_total_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__loss_count',
        'hiloloss_total_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__loss_rt',
        'hiloloss_total_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__loss_std',
        'hiloloss_total_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__loss_count',
        'hiloloss_total_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__loss_rt',
        'hiloloss_total_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__loss_std',
        'hiloloss_pos_feedback_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__loss__crct_count',
        'hiloloss_pos_feedback_bothruns_rate': 'mr_y_tfmri__mid__beh__trial__loss__crct_prop',
        'hiloloss_pos_feedback_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__loss__crct_rt',
        'hiloloss_pos_feedback_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__loss__crct_std',
        'hiloloss_pos_feedback_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__loss__crct_count',
        'hiloloss_pos_feedback_run_1_rate': 'mr_y_tfmri__mid__beh__r01__trial__loss__crct_prop',
        'hiloloss_pos_feedback_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__loss__crct_rt',
        'hiloloss_pos_feedback_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__loss__crct_std',
        'hiloloss_pos_feedback_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__loss__crct_count',
        'hiloloss_pos_feedback_run_2_rate': 'mr_y_tfmri__mid__beh__r02__trial__loss__crct_prop',
        'hiloloss_pos_feedback_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__loss__crct_rt',
        'hiloloss_pos_feedback_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__loss__crct_std',
        'hiloloss_neg_feedback_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__loss__incrct_count',
        'hiloloss_neg_feedback_bothruns_rate': 'mr_y_tfmri__mid__beh__trial__loss__incrct_prop',
        'hiloloss_neg_feedback_bothruns_mean_rt': 'mr_y_tfmri__mid__beh__trial__loss__incrct_rt',
        'hiloloss_neg_feedback_bothruns_std_rt': 'mr_y_tfmri__mid__beh__trial__loss__incrct_std',
        'hiloloss_neg_feedback_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__loss__incrct_count',
        'hiloloss_neg_feedback_run_1_rate': 'mr_y_tfmri__mid__beh__r01__trial__loss__incrct_prop',
        'hiloloss_neg_feedback_run_1_mean_rt': 'mr_y_tfmri__mid__beh__r01__trial__loss__incrct_rt',
        'hiloloss_neg_feedback_run_1_std_rt': 'mr_y_tfmri__mid__beh__r01__trial__loss__incrct_std',
        'hiloloss_neg_feedback_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__loss__incrct_count',
        'hiloloss_neg_feedback_run_2_rate': 'mr_y_tfmri__mid__beh__r02__trial__loss__incrct_prop',
        'hiloloss_neg_feedback_run_2_mean_rt': 'mr_y_tfmri__mid__beh__r02__trial__loss__incrct_rt',
        'hiloloss_neg_feedback_run_2_std_rt': 'mr_y_tfmri__mid__beh__r02__trial__loss__incrct_std',
        'earnings_total': 'mr_y_tfmri__mid__beh__earn_sum',
        'earnings_run_1': 'mr_y_tfmri__mid__beh__r01__earn_sum',
        'earnings_run_2': 'mr_y_tfmri__mid__beh__r02__earn_sum',
        'feedback_flag': 'mr_y_tfmri__mid__beh__feedback_indicator',
        'trialresp_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__resp_count',
        'valid_trialresp_run_1_numtrials': 'mr_y_tfmri__mid__beh__r01__trial__valid__resp_count',
        'trialresp_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__resp_count',
        'valid_trialresp_run_2_numtrials': 'mr_y_tfmri__mid__beh__r02__trial__valid__resp_count',
        'trialresp_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__resp_count',
        'validtrialresp_bothruns_numtrials': 'mr_y_tfmri__mid__beh__trial__valid__resp_count',
    },
    'SST': {
        'switch_flag': 'mr_y_tfmri__sst__beh__switch_indicator',
        'perform_flag': 'mr_y_tfmri__sst__beh__qc_indicator',
        'nruns': 'mr_y_tfmri__sst__beh__run_count',
        'eprime_behav_only': 'mr_y_tfmri__sst__beh__eprime__behav_indicator',
        'eprime_expt_name': 'mr_y_tfmri__sst__beh__eprime__expt_name',
        'eprime_expt_version': 'mr_y_tfmri__sst__beh__eprime__expt_ver',
        'eprime_runtime_version': 'mr_y_tfmri__sst__beh__eprime__runtime_ver',
        'start_delay_min': 'mr_y_tfmri__sst__beh__start__delay_min',
        'start_delay_max': 'mr_y_tfmri__sst__beh__start__delay_max',
        'start_delay_diff_min': 'mr_y_tfmri__sst__beh__start__delay__diff_min',
        'start_delay_diff_max': 'mr_y_tfmri__sst__beh__start__delay__diff_max',
        'init_delay_diff_min': 'mr_y_tfmri__sst__beh__init__delay__diff_min',
        'init_delay_diff_max': 'mr_y_tfmri__sst__beh__init__delay__diff_max',
        'Go_total': 'mr_y_tfmri__sst__beh__go__trial_count',
        'Stop_total': 'mr_y_tfmri__sst__beh__stop__trial_count',
        'Total': 'mr_y_tfmri__sst__beh__trial_count',
        'Go_total_run_1': 'mr_y_tfmri__sst__beh__r01__go__trial_count',
        'Stop_total_run_1': 'mr_y_tfmri__sst__beh__r01__stop__trial_count',
        'Go_total_run_2': 'mr_y_tfmri__sst__beh__r02__go__trial_count',
        'Stop_total_run_2': 'mr_y_tfmri__sst__beh__r02__stop__trial_count',
        'Total_run_1': 'mr_y_tfmri__sst__beh__r01__trial_count',
        'Total_run_2': 'mr_y_tfmri__sst__beh__r02__trial_count',
        'CorrectGo_counts_total': 'mr_y_tfmri__sst__beh__go__crct_count',
        'CorrectGo_counts_run_1': 'mr_y_tfmri__sst__beh__r01__go__crct_count',
        'CorrectGo_counts_run_2': 'mr_y_tfmri__sst__beh__r02__go__crct_count',
        'CorrectGo_percent_total': 'mr_y_tfmri__sst__beh__go__crct_prop',
        'CorrectGo_percent_run_1': 'mr_y_tfmri__sst__beh__r01__go__crct_prop',
        'CorrectGo_percent_run_2': 'mr_y_tfmri__sst__beh__r02__go__crct_prop',
        'CorrectGo_rt_total': 'mr_y_tfmri__sst__beh__go__crct_rt',
        'CorrectGo_rt_std_total': 'mr_y_tfmri__sst__beh__go__crct_std',
        'CorrectGo_rt_run_1': 'mr_y_tfmri__sst__beh__r01__go__crct_rt',
        'CorrectGo_rt_std_run_1': 'mr_y_tfmri__sst__beh__r01__go__crct_std',
        'CorrectGo_rt_run_2': 'mr_y_tfmri__sst__beh__r02__go__crct_rt',
        'CorrectGo_rt_std_run_2': 'mr_y_tfmri__sst__beh__r02__go__crct_std',
        'CorrectLateGo_counts_total': 'mr_y_tfmri__sst__beh__go__crct__late_count',
        'CorrectLateGo_counts_run_1': 'mr_y_tfmri__sst__beh__r01__go__crct__late_count',
        'CorrectLateGo_counts_run_2': 'mr_y_tfmri__sst__beh__r02__go__crct__late_count',
        'CorrectLateGo_percent_total': 'mr_y_tfmri__sst__beh__go__crct__late_prop',
        'CorrectLateGo_percent_run_1': 'mr_y_tfmri__sst__beh__r01__go__crct__late_prop',
        'CorrectLateGo_percent_run_2': 'mr_y_tfmri__sst__beh__r02__go__crct__late_prop',
        'CorrectLateGo_rt_total': 'mr_y_tfmri__sst__beh__go__crct__late_rt',
        'CorrectLateGo_rt_std_total': 'mr_y_tfmri__sst__beh__go__crct__late_std',
        'CorrectLateGo_rt_run_1': 'mr_y_tfmri__sst__beh__r01__go__crct__late_rt',
        'CorrectLateGo_rt_std_run_1': 'mr_y_tfmri__sst__beh__r01__go__crct__late_std',
        'CorrectLateGo_rt_run_2': 'mr_y_tfmri__sst__beh__r02__go__crct__late_rt',
        'CorrectLateGo_rt_std_run_2': 'mr_y_tfmri__sst__beh__r02__go__crct__late_std',
        'IncorrectGo_counts_total': 'mr_y_tfmri__sst__beh__go__incrct_count',
        'IncorrectGo_counts_run_1': 'mr_y_tfmri__sst__beh__r01__go__incrct_count',
        'IncorrectGo_counts_run_2': 'mr_y_tfmri__sst__beh__r02__go__incrct_count',
        'IncorrectGo_percent_total': 'mr_y_tfmri__sst__beh__go__incrct_prop',
        'IncorrectGo_percent_run_1': 'mr_y_tfmri__sst__beh__r01__go__incrct_prop',
        'IncorrectGo_percent_run_2': 'mr_y_tfmri__sst__beh__r02__go__incrct_prop',
        'IncorrectGo_rt_total': 'mr_y_tfmri__sst__beh__go__incrct_rt',
        'IncorrectGo_rt_std_total': 'mr_y_tfmri__sst__beh__go__incrct_std',
        'IncorrectGo_rt_run_1': 'mr_y_tfmri__sst__beh__r01__go__incrct_rt',
        'IncorrectGo_rt_std_run_1': 'mr_y_tfmri__sst__beh__r01__go__incrct_std',
        'IncorrectGo_rt_run_2': 'mr_y_tfmri__sst__beh__r02__go__incrct_rt',
        'IncorrectGo_rt_std_run_2': 'mr_y_tfmri__sst__beh__r02__go__incrct_std',
        'IncorrectLateGo_counts_total': 'mr_y_tfmri__sst__beh__go__incrct__late_count',
        'IncorrectLateGo_counts_run_1': 'mr_y_tfmri__sst__beh__r01__go__incrct__late_count',
        'IncorrectLateGo_counts_run_2': 'mr_y_tfmri__sst__beh__r02__go__incrct__late_count',
        'IncorrectLateGo_percent_total': 'mr_y_tfmri__sst__beh__go__incrct__late_prop',
        'IncorrectLateGo_percent_run_1': 'mr_y_tfmri__sst__beh__r01__go__incrct__late_prop',
        'IncorrectLateGo_percent_run_2': 'mr_y_tfmri__sst__beh__r02__go__incrct__late_prop',
        'IncorrectLateGo_rt_total': 'mr_y_tfmri__sst__beh__go__incrct__late_rt',
        'IncorrectLateGo_rt_std_total': 'mr_y_tfmri__sst__beh__go__incrct__late_std',
        'IncorrectLateGo_rt_run_1': 'mr_y_tfmri__sst__beh__r01__go__incrct__late_rt',
        'IncorrectLateGo_rt_std_run_1': 'mr_y_tfmri__sst__beh__r01__go__incrct__late_std',
        'IncorrectLateGo_rt_run_2': 'mr_y_tfmri__sst__beh__r02__go__incrct__late_rt',
        'IncorrectLateGo_rt_std_run_2': 'mr_y_tfmri__sst__beh__r02__go__incrct__late_std',
        'NoRespGo_counts_total': 'mr_y_tfmri__sst__beh__go__unans_count',
        'NoRespGo_counts_run_1': 'mr_y_tfmri__sst__beh__r01__go__unans_count',
        'NoRespGo_counts_run_2': 'mr_y_tfmri__sst__beh__r02__go__unans_count',
        'NoRespGo_percent_total': 'mr_y_tfmri__sst__beh__go__unans_prop',
        'NoRespGo_percent_run_1': 'mr_y_tfmri__sst__beh__r01__go__unans_prop',
        'NoRespGo_percent_run_2': 'mr_y_tfmri__sst__beh__r02__go__unans_prop',
        'CorrectStop_counts_total': 'mr_y_tfmri__sst__beh__stop__crct_count',
        'CorrectStop_counts_run_1': 'mr_y_tfmri__sst__beh__r01__stop__crct_count',
        'CorrectStop_counts_run_2': 'mr_y_tfmri__sst__beh__r02__stop__crct_count',
        'CorrectStop_percent_total': 'mr_y_tfmri__sst__beh__stop__crct_prop',
        'CorrectStop_percent_run_1': 'mr_y_tfmri__sst__beh__r01__stop__crct_prop',
        'CorrectStop_percent_run_2': 'mr_y_tfmri__sst__beh__r02__stop__crct_prop',
        'IncorrectStop_counts_total': 'mr_y_tfmri__sst__beh__stop__incrct_count',
        'IncorrectStop_counts_run_1': 'mr_y_tfmri__sst__beh__r01__stop__incrct_count',
        'IncorrectStop_counts_run_2': 'mr_y_tfmri__sst__beh__r02__stop__incrct_count',
        'IncorrectStop_percent_total': 'mr_y_tfmri__sst__beh__stop__incrct_prop',
        'IncorrectStop_percent_run_1': 'mr_y_tfmri__sst__beh__r01__stop__incrct_prop',
        'IncorrectStop_percent_run_2': 'mr_y_tfmri__sst__beh__r02__stop__incrct_prop',
        'IncorrectStop_rt_total': 'mr_y_tfmri__sst__beh__stop__incrct_rt',
        'IncorrectStop_rt_std_total': 'mr_y_tfmri__sst__beh__stop__incrct_std',
        'IncorrectStop_rt_run_1': 'mr_y_tfmri__sst__beh__r01__stop__incrct_rt',
        'IncorrectStop_rt_std_run_1': 'mr_y_tfmri__sst__beh__r01__stop__incrct_std',
        'IncorrectStop_rt_run_2': 'mr_y_tfmri__sst__beh__r02__stop__incrct_rt',
        'IncorrectStop_rt_std_run_2': 'mr_y_tfmri__sst__beh__r02__stop__incrct_std',
        'SSD_mean_total': 'mr_y_tfmri__sst__beh__stop__ssd_mean',
        'STE_counts_total': 'mr_y_tfmri__sst__beh__stop__early_count',
        'STE_counts_run_1': 'mr_y_tfmri__sst__beh__r01__stop__early_count',
        'STE_counts_run_2': 'mr_y_tfmri__sst__beh__r02__stop__early_count',
        'STE_percent_total': 'mr_y_tfmri__sst__beh__stop__early_prop',
        'STE_percent_run_1': 'mr_y_tfmri__sst__beh__r01__stop__early_prop',
        'STE_percent_run_2': 'mr_y_tfmri__sst__beh__r02__stop__early_prop',
        'ZeroSSDResp_counts_total': 'mr_y_tfmri__sst__beh__0ssd__resp_count',
        'ZeroSSDResp_counts_run_1': 'mr_y_tfmri__sst__beh__r01__0ssd__resp_count',
        'ZeroSSDResp_counts_run_2': 'mr_y_tfmri__sst__beh__r02__0ssd__resp_count',
        'SSRT_mean_total': 'mr_y_tfmri__sst__beh__ssrt_mean',
        'SSRT_integrated_total': 'mr_y_tfmri__sst__beh__ssrt_intgr',
        'violatorflag': 'mr_y_tfmri__sst__beh__violat_indicator',
        'glitchflag': 'mr_y_tfmri__sst__beh__coderr_indicator',
        'zeroSSDcount': 'mr_y_tfmri__sst__beh__stop__0ssd_count',
    },
    'NBACK': {
        'switch_flag': 'mr_y_tfmri__nback__beh__switch_indicator',
        'perform_flag': 'mr_y_tfmri__nback__beh__qc_indicator',
        'nruns': 'mr_y_tfmri__nback__beh__run_count',
        'eprime_behav_only': 'mr_y_tfmri__nback__beh__eprime__behav_indicator',
        'eprime_expt_name': 'mr_y_tfmri__nback__beh__eprime__expt_name',
        'eprime_expt_version': 'mr_y_tfmri__nback__beh__eprime__expt_ver',
        'eprime_runtime_version': 'mr_y_tfmri__nback__beh__eprime__runtime_ver',
        'start_delay_min': 'mr_y_tfmri__nback__beh__start__delay_min',
        'start_delay_max': 'mr_y_tfmri__nback__beh__start__delay_max',
        'start_delay_diff_min': 'mr_y_tfmri__nback__beh__start__delay__diff_min',
        'start_delay_diff_max': 'mr_y_tfmri__nback__beh__start__delay__diff_max',
        'init_delay_diff_min': 'mr_y_tfmri__nback__beh__init__delay__diff_min',
        'init_delay_diff_max': 'mr_y_tfmri__nback__beh__init__delay__diff_max',
        'all_total': 'mr_y_tfmri__nback__beh__trial_count',
        'all_run_1': 'mr_y_tfmri__nback__beh__r01__trial_count',
        'all_run_2': 'mr_y_tfmri__nback__beh__r02__trial_count',
        'all_total_correct': 'mr_y_tfmri__nback__beh__crct_count',
        'all_total_acc': 'mr_y_tfmri__nback__beh_acc',
        'all_total_correct_rt': 'mr_y_tfmri__nback__beh_rt',
        'all_total_correct_std': 'mr_y_tfmri__nback__beh_std',
        'block_2_back_total': 'mr_y_tfmri__nback__beh__2b_count',
        'block_2_back_run_1': 'mr_y_tfmri__nback__beh__r01__2b_count',
        'block_2_back_run_2': 'mr_y_tfmri__nback__beh__r02__2b_count',
        'block_2_back_total_correct': 'mr_y_tfmri__nback__beh__2b__crct_count',
        'block_2_back_total_acc': 'mr_y_tfmri__nback__beh__2b_acc',
        'block_2_back_total_correct_rt': 'mr_y_tfmri__nback__beh__2b_rt',
        'block_2_back_total_correct_std': 'mr_y_tfmri__nback__beh__2b_std',
        'block_2_back_run_1_correct': 'mr_y_tfmri__nback__beh__r01__2b__crct_count',
        'block_2_back_run_1_acc': 'mr_y_tfmri__nback__beh__r01__2b_acc',
        'block_2_back_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__2b_rt',
        'block_2_back_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__2b_std',
        'block_2_back_run_2_correct': 'mr_y_tfmri__nback__beh__r02__2b__crct_count',
        'block_2_back_run_2_acc': 'mr_y_tfmri__nback__beh__r02__2b_acc',
        'block_2_back_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__2b_rt',
        'block_2_back_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__2b_std',
        'posface_total': 'mr_y_tfmri__nback__beh__pos_count',
        'posface_run_1': 'mr_y_tfmri__nback__beh__r01__pos_count',
        'posface_run_2': 'mr_y_tfmri__nback__beh__r02__pos_count',
        'posface_correct': 'mr_y_tfmri__nback__beh__pos__crct_count',
        'posface_acc': 'mr_y_tfmri__nback__beh__pos_acc',
        'posface_correct_rt': 'mr_y_tfmri__nback__beh__pos_rt',
        'posface_correct_std': 'mr_y_tfmri__nback__beh__pos_std',
        'posface_run_1_correct': 'mr_y_tfmri__nback__beh__r01__pos__crct_count',
        'posface_run_1_acc': 'mr_y_tfmri__nback__beh__r01__pos_acc',
        'posface_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__pos_rt',
        'posface_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__pos_std',
        'posface_run_2_correct': 'mr_y_tfmri__nback__beh__r02__pos__crct_count',
        'posface_run_2_acc': 'mr_y_tfmri__nback__beh__r02__pos_acc',
        'posface_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__pos_rt',
        'posface_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__pos_std',
        'block_2_back_posface_total': 'mr_y_tfmri__nback__beh__2b__pos_count',
        'block_2_back_posface_run_1': 'mr_y_tfmri__nback__beh__r01__2b__pos_count',
        'block_2_back_posface_run_2': 'mr_y_tfmri__nback__beh__r02__2b__pos_count',
        'block_2_back_posface_correct': 'mr_y_tfmri__nback__beh__2b__pos__crct_count',
        'block_2_back_posface_acc': 'mr_y_tfmri__nback__beh__2b__pos_acc',
        'block_2_back_posface_correct_rt': 'mr_y_tfmri__nback__beh__2b__pos_rt',
        'block_2_back_posface_correct_std': 'mr_y_tfmri__nback__beh__2b__pos_std',
        'block_2_back_posface_run_1_correct': 'mr_y_tfmri__nback__beh__r01__2b__pos__crct_count',
        'block_2_back_posface_run_1_acc': 'mr_y_tfmri__nback__beh__r01__2b__pos_acc',
        'block_2_back_posface_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__2b__pos_rt',
        'block_2_back_posface_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__2b__pos_std',
        'block_2_back_posface_run_2_correct': 'mr_y_tfmri__nback__beh__r02__2b__pos__crct_count',
        'block_2_back_posface_run_2_acc': 'mr_y_tfmri__nback__beh__r02__2b__pos_acc',
        'block_2_back_posface_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__2b__pos_rt',
        'block_2_back_posface_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__2b__pos_std',
        'block_2_back_posface_target_total': 'mr_y_tfmri__nback__beh__2b__pos__target_count',
        'block_2_back_posface_target_correct': 'mr_y_tfmri__nback__beh__2b__pos__target__crct_count',
        'block_2_back_posface_target_acc': 'mr_y_tfmri__nback__beh__2b__pos__target_acc',
        'block_2_back_posface_target_rt': 'mr_y_tfmri__nback__beh__2b__pos__target_rt',
        'block_2_back_posface_target_std': 'mr_y_tfmri__nback__beh__2b__pos__target_std',
        'block_2_back_posface_target_run_1_total': 'mr_y_tfmri__nback__beh__r01__2b__pos__target_count',
        'block_2_back_posface_target_run_1_correct': 'mr_y_tfmri__nback__beh__r01__2b__pos__target__crct_count',
        'block_2_back_posface_target_run_1_acc': 'mr_y_tfmri__nback__beh__r01__2b__pos__target_acc',
        'block_2_back_posface_target_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__2b__pos__target_rt',
        'block_2_back_posface_target_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__2b__pos__target_std',
        'block_2_back_posface_target_run_2_total': 'mr_y_tfmri__nback__beh__r02__2b__pos__target_count',
        'block_2_back_posface_target_run_2_correct': 'mr_y_tfmri__nback__beh__r02__2b__pos__target__crct_count',
        'block_2_back_posface_target_run_2_acc': 'mr_y_tfmri__nback__beh__r02__2b__pos__target_acc',
        'block_2_back_posface_target_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__2b__pos__target_rt',
        'block_2_back_posface_target_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__2b__pos__target_std',
        'block_2_back_posface_lure_total': 'mr_y_tfmri__nback__beh__2b__pos__lure_count',
        'block_2_back_posface_lure_correct': 'mr_y_tfmri__nback__beh__2b__pos__lure__crct_count',
        'block_2_back_posface_lure_acc': 'mr_y_tfmri__nback__beh__2b__pos__lure_acc',
        'block_2_back_posface_lure_rt': 'mr_y_tfmri__nback__beh__2b__pos__lure_rt',
        'block_2_back_posface_lure_std': 'mr_y_tfmri__nback__beh__2b__pos__lure_std',
        'block_2_back_posface_lure_run_1_total': 'mr_y_tfmri__nback__beh__r01__2b__pos__lure_count',
        'block_2_back_posface_lure_run_1_correct': 'mr_y_tfmri__nback__beh__r01__2b__pos__lure__crct_count',
        'block_2_back_posface_lure_run_1_acc': 'mr_y_tfmri__nback__beh__r01__2b__pos__lure_acc',
        'block_2_back_posface_lure_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__2b__pos__lure_rt',
        'block_2_back_posface_lure_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__2b__pos__lure_std',
        'block_2_back_posface_lure_run_2_total': 'mr_y_tfmri__nback__beh__r02__2b__pos__lure_count',
        'block_2_back_posface_lure_run_2_correct': 'mr_y_tfmri__nback__beh__r02__2b__pos__lure__crct_count',
        'block_2_back_posface_lure_run_2_acc': 'mr_y_tfmri__nback__beh__r02__2b__pos__lure_acc',
        'block_2_back_posface_lure_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__2b__pos__lure_rt',
        'block_2_back_posface_lure_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__2b__pos__lure_std',
        'block_2_back_posface_nonlure_total': 'mr_y_tfmri__nback__beh__2b__pos__nonlure_count',
        'block_2_back_posface_nonlure_correct': 'mr_y_tfmri__nback__beh__2b__pos__nonlure__crct_count',
        'block_2_back_posface_nonlure_acc': 'mr_y_tfmri__nback__beh__2b__pos__nonlure_acc',
        'block_2_back_posface_nonlure_rt': 'mr_y_tfmri__nback__beh__2b__pos__nonlure_rt',
        'block_2_back_posface_nonlure_std': 'mr_y_tfmri__nback__beh__2b__pos__nonlure_std',
        'block_2_back_posface_nonlure_run_1_total': 'mr_y_tfmri__nback__beh__r01__2b__pos__nonlure_count',
        'block_2_back_posface_nonlure_run_1_correct': 'mr_y_tfmri__nback__beh__r01__2b__pos__nonlure__crct_count',
        'block_2_back_posface_nonlure_run_1_acc': 'mr_y_tfmri__nback__beh__r01__2b__pos__nonlure_acc',
        'block_2_back_posface_nonlure_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__2b__pos__nonlure_rt',
        'block_2_back_posface_nonlure_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__2b__pos__nonlure_std',
        'block_2_back_posface_nonlure_run_2_total': 'mr_y_tfmri__nback__beh__r02__2b__pos__nonlure_count',
        'block_2_back_posface_nonlure_run_2_correct': 'mr_y_tfmri__nback__beh__r02__2b__pos__nonlure__crct_count',
        'block_2_back_posface_nonlure_run_2_acc': 'mr_y_tfmri__nback__beh__r02__2b__pos__nonlure_acc',
        'block_2_back_posface_nonlure_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__2b__pos__nonlure_rt',
        'block_2_back_posface_nonlure_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__2b__pos__nonlure_std',
        'neutface_total': 'mr_y_tfmri__nback__beh__neut_count',
        'neutface_run_1': 'mr_y_tfmri__nback__beh__r01__neut_count',
        'neutface_run_2': 'mr_y_tfmri__nback__beh__r02__neut_count',
        'neutface_correct': 'mr_y_tfmri__nback__beh__neut__crct_count',
        'neutface_acc': 'mr_y_tfmri__nback__beh__neut_acc',
        'neutface_correct_rt': 'mr_y_tfmri__nback__beh__neut_rt',
        'neutface_correct_std': 'mr_y_tfmri__nback__beh__neut_std',
        'neutface_run_1_correct': 'mr_y_tfmri__nback__beh__r01__neut__crct_count',
        'neutface_run_1_acc': 'mr_y_tfmri__nback__beh__r01__neut_acc',
        'neutface_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__neut_rt',
        'neutface_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__neut_std',
        'neutface_run_2_correct': 'mr_y_tfmri__nback__beh__r02__neut__crct_count',
        'neutface_run_2_acc': 'mr_y_tfmri__nback__beh__r02__neut_acc',
        'neutface_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__neut_rt',
        'neutface_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__neut_std',
        'block_2_back_neutface_total': 'mr_y_tfmri__nback__beh__2b__neut_count',
        'block_2_back_neutface_run_1': 'mr_y_tfmri__nback__beh__r01__2b__neut_count',
        'block_2_back_neutface_run_2': 'mr_y_tfmri__nback__beh__r02__2b__neut_count',
        'block_2_back_neutface_correct': 'mr_y_tfmri__nback__beh__2b__neut__crct_count',
        'block_2_back_neutface_acc': 'mr_y_tfmri__nback__beh__2b__neut_acc',
        'block_2_back_neutface_correct_rt': 'mr_y_tfmri__nback__beh__2b__neut_rt',
        'block_2_back_neutface_correct_std': 'mr_y_tfmri__nback__beh__2b__neut_std',
        'block_2_back_neutface_run_1_correct': 'mr_y_tfmri__nback__beh__r01__2b__neut__crct_count',
        'block_2_back_neutface_run_1_acc': 'mr_y_tfmri__nback__beh__r01__2b__neut_acc',
        'block_2_back_neutface_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__2b__neut_rt',
        'block_2_back_neutface_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__2b__neut_std',
        'block_2_back_neutface_run_2_correct': 'mr_y_tfmri__nback__beh__r02__2b__neut__crct_count',
        'block_2_back_neutface_run_2_acc': 'mr_y_tfmri__nback__beh__r02__2b__neut_acc',
        'block_2_back_neutface_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__2b__neut_rt',
        'block_2_back_neutface_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__2b__neut_std',
        'block_2_back_neutface_target_total': 'mr_y_tfmri__nback__beh__2b__neut__target_count',
        'block_2_back_neutface_target_correct': 'mr_y_tfmri__nback__beh__2b__neut__target__crct_count',
        'block_2_back_neutface_target_acc': 'mr_y_tfmri__nback__beh__2b__neut__target_acc',
        'block_2_back_neutface_target_rt': 'mr_y_tfmri__nback__beh__2b__neut__target_rt',
        'block_2_back_neutface_target_std': 'mr_y_tfmri__nback__beh__2b__neut__target_std',
        'block_2_back_neutface_target_run_1_total': 'mr_y_tfmri__nback__beh__r01__2b__neut__target_count',
        'block_2_back_neutface_target_run_1_correct': 'mr_y_tfmri__nback__beh__r01__2b__neut__target__crct_count',
        'block_2_back_neutface_target_run_1_acc': 'mr_y_tfmri__nback__beh__r01__2b__neut__target_acc',
        'block_2_back_neutface_target_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__2b__neut__target_rt',
        'block_2_back_neutface_target_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__2b__neut__target_std',
        'block_2_back_neutface_target_run_2_total': 'mr_y_tfmri__nback__beh__r02__2b__neut__target_count',
        'block_2_back_neutface_target_run_2_correct': 'mr_y_tfmri__nback__beh__r02__2b__neut__target__crct_count',
        'block_2_back_neutface_target_run_2_acc': 'mr_y_tfmri__nback__beh__r02__2b__neut__target_acc',
        'block_2_back_neutface_target_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__2b__neut__target_rt',
        'block_2_back_neutface_target_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__2b__neut__target_std',
        'block_2_back_neutface_lure_total': 'mr_y_tfmri__nback__beh__2b__neut__lure_count',
        'block_2_back_neutface_lure_correct': 'mr_y_tfmri__nback__beh__2b__neut__lure__crct_count',
        'block_2_back_neutface_lure_acc': 'mr_y_tfmri__nback__beh__2b__neut__lure_acc',
        'block_2_back_neutface_lure_rt': 'mr_y_tfmri__nback__beh__2b__neut__lure_rt',
        'block_2_back_neutface_lure_std': 'mr_y_tfmri__nback__beh__2b__neut__lure_std',
        'block_2_back_neutface_lure_run_1_total': 'mr_y_tfmri__nback__beh__r01__2b__neut__lure_count',
        'block_2_back_neutface_lure_run_1_correct': 'mr_y_tfmri__nback__beh__r01__2b__neut__lure__crct_count',
        'block_2_back_neutface_lure_run_1_acc': 'mr_y_tfmri__nback__beh__r01__2b__neut__lure_acc',
        'block_2_back_neutface_lure_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__2b__neut__lure_rt',
        'block_2_back_neutface_lure_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__2b__neut__lure_std',
        'block_2_back_neutface_lure_run_2_total': 'mr_y_tfmri__nback__beh__r02__2b__neut__lure_count',
        'block_2_back_neutface_lure_run_2_correct': 'mr_y_tfmri__nback__beh__r02__2b__neut__lure__crct_count',
        'block_2_back_neutface_lure_run_2_acc': 'mr_y_tfmri__nback__beh__r02__2b__neut__lure_acc',
        'block_2_back_neutface_lure_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__2b__neut__lure_rt',
        'block_2_back_neutface_lure_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__2b__neut__lure_std',
        'block_2_back_neutface_nonlure_total': 'mr_y_tfmri__nback__beh__2b__neut__nonlure_count',
        'block_2_back_neutface_nonlure_correct': 'mr_y_tfmri__nback__beh__2b__neut__nonlure__crct_count',
        'block_2_back_neutface_nonlure_acc': 'mr_y_tfmri__nback__beh__2b__neut__nonlure_acc',
        'block_2_back_neutface_nonlure_rt': 'mr_y_tfmri__nback__beh__2b__neut__nonlure_rt',
        'block_2_back_neutface_nonlure_std': 'mr_y_tfmri__nback__beh__2b__neut__nonlure_std',
        'block_2_back_neutface_nonlure_run_1_total': 'mr_y_tfmri__nback__beh__r01__2b__neut__nonlure_count',
        'block_2_back_neutface_nonlure_run_1_correct': 'mr_y_tfmri__nback__beh__r01__2b__neut__nonlure__crct_count',
        'block_2_back_neutface_nonlure_run_1_acc': 'mr_y_tfmri__nback__beh__r01__2b__neut__nonlure_acc',
        'block_2_back_neutface_nonlure_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__2b__neut__nonlure_rt',
        'block_2_back_neutface_nonlure_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__2b__neut__nonlure_std',
        'block_2_back_neutface_nonlure_run_2_total': 'mr_y_tfmri__nback__beh__r02__2b__neut__nonlure_count',
        'block_2_back_neutface_nonlure_run_2_correct': 'mr_y_tfmri__nback__beh__r02__2b__neut__nonlure__crct_count',
        'block_2_back_neutface_nonlure_run_2_acc': 'mr_y_tfmri__nback__beh__r02__2b__neut__nonlure_acc',
        'block_2_back_neutface_nonlure_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__2b__neut__nonlure_rt',
        'block_2_back_neutface_nonlure_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__2b__neut__nonlure_std',
        'negface_total': 'mr_y_tfmri__nback__beh__neg_count',
        'negface_run_1': 'mr_y_tfmri__nback__beh__r01__neg_count',
        'negface_run_2': 'mr_y_tfmri__nback__beh__r02__neg_count',
        'negface_correct': 'mr_y_tfmri__nback__beh__neg__crct_count',
        'negface_acc': 'mr_y_tfmri__nback__beh__neg_acc',
        'negface_correct_rt': 'mr_y_tfmri__nback__beh__neg_rt',
        'negface_correct_std': 'mr_y_tfmri__nback__beh__neg_std',
        'negface_run_1_correct': 'mr_y_tfmri__nback__beh__r01__neg__crct_count',
        'negface_run_1_acc': 'mr_y_tfmri__nback__beh__r01__neg_acc',
        'negface_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__neg_rt',
        'negface_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__neg_std',
        'negface_run_2_correct': 'mr_y_tfmri__nback__beh__r02__neg__crct_count',
        'negface_run_2_acc': 'mr_y_tfmri__nback__beh__r02__neg_acc',
        'negface_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__neg_rt',
        'negface_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__neg_std',
        'block_2_back_negface_total': 'mr_y_tfmri__nback__beh__2b__neg_count',
        'block_2_back_negface_run_1': 'mr_y_tfmri__nback__beh__r01__2b__neg_count',
        'block_2_back_negface_run_2': 'mr_y_tfmri__nback__beh__r02__2b__neg_count',
        'block_2_back_negface_correct': 'mr_y_tfmri__nback__beh__2b__neg__crct_count',
        'block_2_back_negface_acc': 'mr_y_tfmri__nback__beh__2b__neg_acc',
        'block_2_back_negface_correct_rt': 'mr_y_tfmri__nback__beh__2b__neg_rt',
        'block_2_back_negface_correct_std': 'mr_y_tfmri__nback__beh__2b__neg_std',
        'block_2_back_negface_run_1_correct': 'mr_y_tfmri__nback__beh__r01__2b__neg__crct_count',
        'block_2_back_negface_run_1_acc': 'mr_y_tfmri__nback__beh__r01__2b__neg_acc',
        'block_2_back_negface_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__2b__neg_rt',
        'block_2_back_negface_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__2b__neg_std',
        'block_2_back_negface_run_2_correct': 'mr_y_tfmri__nback__beh__r02__2b__neg__crct_count',
        'block_2_back_negface_run_2_acc': 'mr_y_tfmri__nback__beh__r02__2b__neg_acc',
        'block_2_back_negface_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__2b__neg_rt',
        'block_2_back_negface_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__2b__neg_std',
        'block_2_back_negface_target_total': 'mr_y_tfmri__nback__beh__2b__neg__target_count',
        'block_2_back_negface_target_correct': 'mr_y_tfmri__nback__beh__2b__neg__target__crct_count',
        'block_2_back_negface_target_acc': 'mr_y_tfmri__nback__beh__2b__neg__target_acc',
        'block_2_back_negface_target_rt': 'mr_y_tfmri__nback__beh__2b__neg__target_rt',
        'block_2_back_negface_target_std': 'mr_y_tfmri__nback__beh__2b__neg__target_std',
        'block_2_back_negface_target_run_1_total': 'mr_y_tfmri__nback__beh__r01__2b__neg__target_count',
        'block_2_back_negface_target_run_1_correct': 'mr_y_tfmri__nback__beh__r01__2b__neg__target__crct_count',
        'block_2_back_negface_target_run_1_acc': 'mr_y_tfmri__nback__beh__r01__2b__neg__target_acc',
        'block_2_back_negface_target_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__2b__neg__target_rt',
        'block_2_back_negface_target_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__2b__neg__target_std',
        'block_2_back_negface_target_run_2_total': 'mr_y_tfmri__nback__beh__r02__2b__neg__target_count',
        'block_2_back_negface_target_run_2_correct': 'mr_y_tfmri__nback__beh__r02__2b__neg__target__crct_count',
        'block_2_back_negface_target_run_2_acc': 'mr_y_tfmri__nback__beh__r02__2b__neg__target_acc',
        'block_2_back_negface_target_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__2b__neg__target_rt',
        'block_2_back_negface_target_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__2b__neg__target_std',
        'block_2_back_negface_lure_total': 'mr_y_tfmri__nback__beh__2b__neg__lure_count',
        'block_2_back_negface_lure_correct': 'mr_y_tfmri__nback__beh__2b__neg__lure__crct_count',
        'block_2_back_negface_lure_acc': 'mr_y_tfmri__nback__beh__2b__neg__lure_acc',
        'block_2_back_negface_lure_rt': 'mr_y_tfmri__nback__beh__2b__neg__lure_rt',
        'block_2_back_negface_lure_std': 'mr_y_tfmri__nback__beh__2b__neg__lure_std',
        'block_2_back_negface_lure_run_1_total': 'mr_y_tfmri__nback__beh__r01__2b__neg__lure_count',
        'block_2_back_negface_lure_run_1_correct': 'mr_y_tfmri__nback__beh__r01__2b__neg__lure__crct_count',
        'block_2_back_negface_lure_run_1_acc': 'mr_y_tfmri__nback__beh__r01__2b__neg__lure_acc',
        'block_2_back_negface_lure_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__2b__neg__lure_rt',
        'block_2_back_negface_lure_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__2b__neg__lure_std',
        'block_2_back_negface_lure_run_2_total': 'mr_y_tfmri__nback__beh__r02__2b__neg__lure_count',
        'block_2_back_negface_lure_run_2_correct': 'mr_y_tfmri__nback__beh__r02__2b__neg__lure__crct_count',
        'block_2_back_negface_lure_run_2_acc': 'mr_y_tfmri__nback__beh__r02__2b__neg__lure_acc',
        'block_2_back_negface_lure_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__2b__neg__lure_rt',
        'block_2_back_negface_lure_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__2b__neg__lure_std',
        'block_2_back_negface_nonlure_total': 'mr_y_tfmri__nback__beh__2b__neg__nonlure_count',
        'block_2_back_negface_nonlure_correct': 'mr_y_tfmri__nback__beh__2b__neg__nonlure__crct_count',
        'block_2_back_negface_nonlure_acc': 'mr_y_tfmri__nback__beh__2b__neg__nonlure_acc',
        'block_2_back_negface_nonlure_rt': 'mr_y_tfmri__nback__beh__2b__neg__nonlure_rt',
        'block_2_back_negface_nonlure_std': 'mr_y_tfmri__nback__beh__2b__neg__nonlure_std',
        'block_2_back_negface_nonlure_run_1_total': 'mr_y_tfmri__nback__beh__r01__2b__neg__nonlure_count',
        'block_2_back_negface_nonlure_run_1_correct': 'mr_y_tfmri__nback__beh__r01__2b__neg__nonlure__crct_count',
        'block_2_back_negface_nonlure_run_1_acc': 'mr_y_tfmri__nback__beh__r01__2b__neg__nonlure_acc',
        'block_2_back_negface_nonlure_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__2b__neg__nonlure_rt',
        'block_2_back_negface_nonlure_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__2b__neg__nonlure_std',
        'block_2_back_negface_nonlure_run_2_total': 'mr_y_tfmri__nback__beh__r02__2b__neg__nonlure_count',
        'block_2_back_negface_nonlure_run_2_correct': 'mr_y_tfmri__nback__beh__r02__2b__neg__nonlure__crct_count',
        'block_2_back_negface_nonlure_run_2_acc': 'mr_y_tfmri__nback__beh__r02__2b__neg__nonlure_acc',
        'block_2_back_negface_nonlure_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__2b__neg__nonlure_rt',
        'block_2_back_negface_nonlure_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__2b__neg__nonlure_std',
        'place_total': 'mr_y_tfmri__nback__beh__place_count',
        'place_run_1': 'mr_y_tfmri__nback__beh__r01__place_count',
        'place_run_2': 'mr_y_tfmri__nback__beh__r02__place_count',
        'place_correct': 'mr_y_tfmri__nback__beh__place__crct_count',
        'place_acc': 'mr_y_tfmri__nback__beh__place_acc',
        'place_correct_rt': 'mr_y_tfmri__nback__beh__place_rt',
        'place_correct_std': 'mr_y_tfmri__nback__beh__place_std',
        'place_run_1_correct': 'mr_y_tfmri__nback__beh__r01__place__crct_count',
        'place_run_1_acc': 'mr_y_tfmri__nback__beh__r01__place_acc',
        'place_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__place_rt',
        'place_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__place_std',
        'place_run_2_correct': 'mr_y_tfmri__nback__beh__r02__place__crct_count',
        'place_run_2_acc': 'mr_y_tfmri__nback__beh__r02__place_acc',
        'place_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__place_rt',
        'place_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__place_std',
        'block_2_back_place_total': 'mr_y_tfmri__nback__beh__2b__place_count',
        'block_2_back_place_run_1': 'mr_y_tfmri__nback__beh__r01__2b__place_count',
        'block_2_back_place_run_2': 'mr_y_tfmri__nback__beh__r02__2b__place_count',
        'block_2_back_place_correct': 'mr_y_tfmri__nback__beh__2b__place__crct_count',
        'block_2_back_place_acc': 'mr_y_tfmri__nback__beh__2b__place_acc',
        'block_2_back_place_correct_rt': 'mr_y_tfmri__nback__beh__2b__place_rt',
        'block_2_back_place_correct_std': 'mr_y_tfmri__nback__beh__2b__place_std',
        'block_2_back_place_run_1_correct': 'mr_y_tfmri__nback__beh__r01__2b__place__crct_count',
        'block_2_back_place_run_1_acc': 'mr_y_tfmri__nback__beh__r01__2b__place_acc',
        'block_2_back_place_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__2b__place_rt',
        'block_2_back_place_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__2b__place_std',
        'block_2_back_place_run_2_correct': 'mr_y_tfmri__nback__beh__r02__2b__place__crct_count',
        'block_2_back_place_run_2_acc': 'mr_y_tfmri__nback__beh__r02__2b__place_acc',
        'block_2_back_place_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__2b__place_rt',
        'block_2_back_place_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__2b__place_std',
        'block_2_back_place_target_total': 'mr_y_tfmri__nback__beh__2b__place__target_count',
        'block_2_back_place_target_correct': 'mr_y_tfmri__nback__beh__2b__place__target__crct_count',
        'block_2_back_place_target_acc': 'mr_y_tfmri__nback__beh__2b__place__target_acc',
        'block_2_back_place_target_rt': 'mr_y_tfmri__nback__beh__2b__place__target_rt',
        'block_2_back_place_target_std': 'mr_y_tfmri__nback__beh__2b__place__target_std',
        'block_2_back_place_target_run_1_total': 'mr_y_tfmri__nback__beh__r01__2b__place__target_count',
        'block_2_back_place_target_run_1_correct': 'mr_y_tfmri__nback__beh__r01__2b__place__target__crct_count',
        'block_2_back_place_target_run_1_acc': 'mr_y_tfmri__nback__beh__r01__2b__place__target_acc',
        'block_2_back_place_target_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__2b__place__target_rt',
        'block_2_back_place_target_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__2b__place__target_std',
        'block_2_back_place_target_run_2_total': 'mr_y_tfmri__nback__beh__r02__2b__place__target_count',
        'block_2_back_place_target_run_2_correct': 'mr_y_tfmri__nback__beh__r02__2b__place__target__crct_count',
        'block_2_back_place_target_run_2_acc': 'mr_y_tfmri__nback__beh__r02__2b__place__target_acc',
        'block_2_back_place_target_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__2b__place__target_rt',
        'block_2_back_place_target_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__2b__place__target_std',
        'block_2_back_place_lure_total': 'mr_y_tfmri__nback__beh__2b__place__lure_count',
        'block_2_back_place_lure_correct': 'mr_y_tfmri__nback__beh__2b__place__lure__crct_count',
        'block_2_back_place_lure_acc': 'mr_y_tfmri__nback__beh__2b__place__lure_acc',
        'block_2_back_place_lure_rt': 'mr_y_tfmri__nback__beh__2b__place__lure_rt',
        'block_2_back_place_lure_std': 'mr_y_tfmri__nback__beh__2b__place__lure_std',
        'block_2_back_place_lure_run_1_total': 'mr_y_tfmri__nback__beh__r01__2b__place__lure_count',
        'block_2_back_place_lure_run_1_correct': 'mr_y_tfmri__nback__beh__r01__2b__place__lure__crct_count',
        'block_2_back_place_lure_run_1_acc': 'mr_y_tfmri__nback__beh__r01__2b__place__lure_acc',
        'block_2_back_place_lure_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__2b__place__lure_rt',
        'block_2_back_place_lure_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__2b__place__lure_std',
        'block_2_back_place_lure_run_2_total': 'mr_y_tfmri__nback__beh__r02__2b__place__lure_count',
        'block_2_back_place_lure_run_2_correct': 'mr_y_tfmri__nback__beh__r02__2b__place__lure__crct_count',
        'block_2_back_place_lure_run_2_acc': 'mr_y_tfmri__nback__beh__r02__2b__place__lure_acc',
        'block_2_back_place_lure_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__2b__place__lure_rt',
        'block_2_back_place_lure_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__2b__place__lure_std',
        'block_2_back_place_nonlure_total': 'mr_y_tfmri__nback__beh__2b__place__nonlure_count',
        'block_2_back_place_nonlure_correct': 'mr_y_tfmri__nback__beh__2b__place__nonlure__crct_count',
        'block_2_back_place_nonlure_acc': 'mr_y_tfmri__nback__beh__2b__place__nonlure_acc',
        'block_2_back_place_nonlure_rt': 'mr_y_tfmri__nback__beh__2b__place__nonlure_rt',
        'block_2_back_place_nonlure_std': 'mr_y_tfmri__nback__beh__2b__place__nonlure_std',
        'block_2_back_place_nonlure_run_1_total': 'mr_y_tfmri__nback__beh__r01__2b__place__nonlure_count',
        'block_2_back_place_nonlure_run_1_correct': 'mr_y_tfmri__nback__beh__r01__2b__place__nonlure__crct_count',
        'block_2_back_place_nonlure_run_1_acc': 'mr_y_tfmri__nback__beh__r01__2b__place__nonlure_acc',
        'block_2_back_place_nonlure_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__2b__place__nonlure_rt',
        'block_2_back_place_nonlure_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__2b__place__nonlure_std',
        'block_2_back_place_nonlure_run_2_total': 'mr_y_tfmri__nback__beh__r02__2b__place__nonlure_count',
        'block_2_back_place_nonlure_run_2_correct': 'mr_y_tfmri__nback__beh__r02__2b__place__nonlure__crct_count',
        'block_2_back_place_nonlure_run_2_acc': 'mr_y_tfmri__nback__beh__r02__2b__place__nonlure_acc',
        'block_2_back_place_nonlure_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__2b__place__nonlure_rt',
        'block_2_back_place_nonlure_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__2b__place__nonlure_std',
        'block_0_back_total': 'mr_y_tfmri__nback__beh__0b_count',
        'block_0_back_run_1': 'mr_y_tfmri__nback__beh__r01__0b_count',
        'block_0_back_run_2': 'mr_y_tfmri__nback__beh__r02__0b_count',
        'block_0_back_total_correct': 'mr_y_tfmri__nback__beh__0b__crct_count',
        'block_0_back_total_acc': 'mr_y_tfmri__nback__beh__0b_acc',
        'block_0_back_total_correct_rt': 'mr_y_tfmri__nback__beh__0b_rt',
        'block_0_back_total_correct_std': 'mr_y_tfmri__nback__beh__0b_std',
        'block_0_back_run_1_correct': 'mr_y_tfmri__nback__beh__r01__0b__crct_count',
        'block_0_back_run_1_acc': 'mr_y_tfmri__nback__beh__r01__0b_acc',
        'block_0_back_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__0b_rt',
        'block_0_back_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__0b_std',
        'block_0_back_run_2_correct': 'mr_y_tfmri__nback__beh__r02__0b__crct_count',
        'block_0_back_run_2_acc': 'mr_y_tfmri__nback__beh__r02__0b_acc',
        'block_0_back_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__0b_rt',
        'block_0_back_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__0b_std',
        'block_0_back_posface_total': 'mr_y_tfmri__nback__beh__0b__pos_count',
        'block_0_back_posface_run_1': 'mr_y_tfmri__nback__beh__r01__0b__pos_count',
        'block_0_back_posface_run_2': 'mr_y_tfmri__nback__beh__r02__0b__pos_count',
        'block_0_back_posface_correct': 'mr_y_tfmri__nback__beh__0b__pos__crct_count',
        'block_0_back_posface_acc': 'mr_y_tfmri__nback__beh__0b__pos_acc',
        'block_0_back_posface_correct_rt': 'mr_y_tfmri__nback__beh__0b__pos_rt',
        'block_0_back_posface_correct_std': 'mr_y_tfmri__nback__beh__0b__pos_std',
        'block_0_back_posface_run_1_correct': 'mr_y_tfmri__nback__beh__r01__0b__pos__crct_count',
        'block_0_back_posface_run_1_acc': 'mr_y_tfmri__nback__beh__r01__0b__pos_acc',
        'block_0_back_posface_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__0b__pos_rt',
        'block_0_back_posface_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__0b__pos_std',
        'block_0_back_posface_run_2_correct': 'mr_y_tfmri__nback__beh__r02__0b__pos__crct_count',
        'block_0_back_posface_run_2_acc': 'mr_y_tfmri__nback__beh__r02__0b__pos_acc',
        'block_0_back_posface_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__0b__pos_rt',
        'block_0_back_posface_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__0b__pos_std',
        'block_0_back_posface_target_total': 'mr_y_tfmri__nback__beh__0b__pos__target_count',
        'block_0_back_posface_target_correct': 'mr_y_tfmri__nback__beh__0b__pos__target__crct_count',
        'block_0_back_posface_target_acc': 'mr_y_tfmri__nback__beh__0b__pos__target_acc',
        'block_0_back_posface_target_rt': 'mr_y_tfmri__nback__beh__0b__pos__target_rt',
        'block_0_back_posface_target_std': 'mr_y_tfmri__nback__beh__0b__pos__target_std',
        'block_0_back_posface_target_run_1_total': 'mr_y_tfmri__nback__beh__r01__0b__pos__target_count',
        'block_0_back_posface_target_run_1_correct': 'mr_y_tfmri__nback__beh__r01__0b__pos__target__crct_count',
        'block_0_back_posface_target_run_1_acc': 'mr_y_tfmri__nback__beh__r01__0b__pos__target_acc',
        'block_0_back_posface_target_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__0b__pos__target_rt',
        'block_0_back_posface_target_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__0b__pos__target_std',
        'block_0_back_posface_target_run_2_total': 'mr_y_tfmri__nback__beh__r02__0b__pos__target_count',
        'block_0_back_posface_target_run_2_correct': 'mr_y_tfmri__nback__beh__r02__0b__pos__target__crct_count',
        'block_0_back_posface_target_run_2_acc': 'mr_y_tfmri__nback__beh__r02__0b__pos__target_acc',
        'block_0_back_posface_target_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__0b__pos__target_rt',
        'block_0_back_posface_target_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__0b__pos__target_std',
        'block_0_back_posface_lure_total': 'mr_y_tfmri__nback__beh__0b__pos__lure_count',
        'block_0_back_posface_lure_correct': 'mr_y_tfmri__nback__beh__0b__pos__lure__crct_count',
        'block_0_back_posface_lure_acc': 'mr_y_tfmri__nback__beh__0b__pos__lure_acc',
        'block_0_back_posface_lure_rt': 'mr_y_tfmri__nback__beh__0b__pos__lure_rt',
        'block_0_back_posface_lure_std': 'mr_y_tfmri__nback__beh__0b__pos__lure_std',
        'block_0_back_posface_lure_run_1_total': 'mr_y_tfmri__nback__beh__r01__0b__pos__lure_count',
        'block_0_back_posface_lure_run_1_correct': 'mr_y_tfmri__nback__beh__r01__0b__pos__lure__crct_count',
        'block_0_back_posface_lure_run_1_acc': 'mr_y_tfmri__nback__beh__r01__0b__pos__lure_acc',
        'block_0_back_posface_lure_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__0b__pos__lure_rt',
        'block_0_back_posface_lure_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__0b__pos__lure_std',
        'block_0_back_posface_lure_run_2_total': 'mr_y_tfmri__nback__beh__r02__0b__pos__lure_count',
        'block_0_back_posface_lure_run_2_correct': 'mr_y_tfmri__nback__beh__r02__0b__pos__lure__crct_count',
        'block_0_back_posface_lure_run_2_acc': 'mr_y_tfmri__nback__beh__r02__0b__pos__lure_acc',
        'block_0_back_posface_lure_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__0b__pos__lure_rt',
        'block_0_back_posface_lure_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__0b__pos__lure_std',
        'block_0_back_posface_nonlure_total': 'mr_y_tfmri__nback__beh__0b__pos__nonlure_count',
        'block_0_back_posface_nonlure_correct': 'mr_y_tfmri__nback__beh__0b__pos__nonlure__crct_count',
        'block_0_back_posface_nonlure_acc': 'mr_y_tfmri__nback__beh__0b__pos__nonlure_acc',
        'block_0_back_posface_nonlure_rt': 'mr_y_tfmri__nback__beh__0b__pos__nonlure_rt',
        'block_0_back_posface_nonlure_std': 'mr_y_tfmri__nback__beh__0b__pos__nonlure_std',
        'block_0_back_posface_nonlure_run_1_total': 'mr_y_tfmri__nback__beh__r01__0b__pos__nonlure_count',
        'block_0_back_posface_nonlure_run_1_correct': 'mr_y_tfmri__nback__beh__r01__0b__pos__nonlure__crct_count',
        'block_0_back_posface_nonlure_run_1_acc': 'mr_y_tfmri__nback__beh__r01__0b__pos__nonlure_acc',
        'block_0_back_posface_nonlure_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__0b__pos__nonlure_rt',
        'block_0_back_posface_nonlure_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__0b__pos__nonlure_std',
        'block_0_back_posface_nonlure_run_2_total': 'mr_y_tfmri__nback__beh__r02__0b__pos__nonlure_count',
        'block_0_back_posface_nonlure_run_2_correct': 'mr_y_tfmri__nback__beh__r02__0b__pos__nonlure__crct_count',
        'block_0_back_posface_nonlure_run_2_acc': 'mr_y_tfmri__nback__beh__r02__0b__pos__nonlure_acc',
        'block_0_back_posface_nonlure_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__0b__pos__nonlure_rt',
        'block_0_back_posface_nonlure_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__0b__pos__nonlure_std',
        'block_0_back_neutface_total': 'mr_y_tfmri__nback__beh__0b__neut_count',
        'block_0_back_neutface_run_1': 'mr_y_tfmri__nback__beh__r01__0b__neut_count',
        'block_0_back_neutface_run_2': 'mr_y_tfmri__nback__beh__r02__0b__neut_count',
        'block_0_back_neutface_correct': 'mr_y_tfmri__nback__beh__0b__neut__crct_count',
        'block_0_back_neutface_acc': 'mr_y_tfmri__nback__beh__0b__neut_acc',
        'block_0_back_neutface_correct_rt': 'mr_y_tfmri__nback__beh__0b__neut_rt',
        'block_0_back_neutface_correct_std': 'mr_y_tfmri__nback__beh__0b__neut_std',
        'block_0_back_neutface_run_1_correct': 'mr_y_tfmri__nback__beh__r01__0b__neut__crct_count',
        'block_0_back_neutface_run_1_acc': 'mr_y_tfmri__nback__beh__r01__0b__neut_acc',
        'block_0_back_neutface_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__0b__neut_rt',
        'block_0_back_neutface_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__0b__neut_std',
        'block_0_back_neutface_run_2_correct': 'mr_y_tfmri__nback__beh__r02__0b__neut__crct_count',
        'block_0_back_neutface_run_2_acc': 'mr_y_tfmri__nback__beh__r02__0b__neut_acc',
        'block_0_back_neutface_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__0b__neut_rt',
        'block_0_back_neutface_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__0b__neut_std',
        'block_0_back_neutface_target_total': 'mr_y_tfmri__nback__beh__0b__neut__target_count',
        'block_0_back_neutface_target_correct': 'mr_y_tfmri__nback__beh__0b__neut__target__crct_count',
        'block_0_back_neutface_target_acc': 'mr_y_tfmri__nback__beh__0b__neut__target_acc',
        'block_0_back_neutface_target_rt': 'mr_y_tfmri__nback__beh__0b__neut__target_rt',
        'block_0_back_neutface_target_std': 'mr_y_tfmri__nback__beh__0b__neut__target_std',
        'block_0_back_neutface_target_run_1_total': 'mr_y_tfmri__nback__beh__r01__0b__neut__target_count',
        'block_0_back_neutface_target_run_1_correct': 'mr_y_tfmri__nback__beh__r01__0b__neut__target__crct_count',
        'block_0_back_neutface_target_run_1_acc': 'mr_y_tfmri__nback__beh__r01__0b__neut__target_acc',
        'block_0_back_neutface_target_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__0b__neut__target_rt',
        'block_0_back_neutface_target_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__0b__neut__target_std',
        'block_0_back_neutface_target_run_2_total': 'mr_y_tfmri__nback__beh__r02__0b__neut__target_count',
        'block_0_back_neutface_target_run_2_correct': 'mr_y_tfmri__nback__beh__r02__0b__neut__target__crct_count',
        'block_0_back_neutface_target_run_2_acc': 'mr_y_tfmri__nback__beh__r02__0b__neut__target_acc',
        'block_0_back_neutface_target_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__0b__neut__target_rt',
        'block_0_back_neutface_target_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__0b__neut__target_std',
        'block_0_back_neutface_lure_total': 'mr_y_tfmri__nback__beh__0b__neut__lure_count',
        'block_0_back_neutface_lure_correct': 'mr_y_tfmri__nback__beh__0b__neut__lure__crct_count',
        'block_0_back_neutface_lure_acc': 'mr_y_tfmri__nback__beh__0b__neut__lure_acc',
        'block_0_back_neutface_lure_rt': 'mr_y_tfmri__nback__beh__0b__neut__lure_rt',
        'block_0_back_neutface_lure_std': 'mr_y_tfmri__nback__beh__0b__neut__lure_std',
        'block_0_back_neutface_lure_run_1_total': 'mr_y_tfmri__nback__beh__r01__0b__neut__lure_count',
        'block_0_back_neutface_lure_run_1_correct': 'mr_y_tfmri__nback__beh__r01__0b__neut__lure__crct_count',
        'block_0_back_neutface_lure_run_1_acc': 'mr_y_tfmri__nback__beh__r01__0b__neut__lure_acc',
        'block_0_back_neutface_lure_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__0b__neut__lure_rt',
        'block_0_back_neutface_lure_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__0b__neut__lure_std',
        'block_0_back_neutface_lure_run_2_total': 'mr_y_tfmri__nback__beh__r02__0b__neut__lure_count',
        'block_0_back_neutface_lure_run_2_correct': 'mr_y_tfmri__nback__beh__r02__0b__neut__lure__crct_count',
        'block_0_back_neutface_lure_run_2_acc': 'mr_y_tfmri__nback__beh__r02__0b__neut__lure_acc',
        'block_0_back_neutface_lure_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__0b__neut__lure_rt',
        'block_0_back_neutface_lure_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__0b__neut__lure_std',
        'block_0_back_neutface_nonlure_total': 'mr_y_tfmri__nback__beh__0b__neut__nonlure_count',
        'block_0_back_neutface_nonlure_correct': 'mr_y_tfmri__nback__beh__0b__neut__nonlure__crct_count',
        'block_0_back_neutface_nonlure_acc': 'mr_y_tfmri__nback__beh__0b__neut__nonlure_acc',
        'block_0_back_neutface_nonlure_rt': 'mr_y_tfmri__nback__beh__0b__neut__nonlure_rt',
        'block_0_back_neutface_nonlure_std': 'mr_y_tfmri__nback__beh__0b__neut__nonlure_std',
        'block_0_back_neutface_nonlure_run_1_total': 'mr_y_tfmri__nback__beh__r01__0b__neut__nonlure_count',
        'block_0_back_neutface_nonlure_run_1_correct': 'mr_y_tfmri__nback__beh__r01__0b__neut__nonlure__crct_count',
        'block_0_back_neutface_nonlure_run_1_acc': 'mr_y_tfmri__nback__beh__r01__0b__neut__nonlure_acc',
        'block_0_back_neutface_nonlure_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__0b__neut__nonlure_rt',
        'block_0_back_neutface_nonlure_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__0b__neut__nonlure_std',
        'block_0_back_neutface_nonlure_run_2_total': 'mr_y_tfmri__nback__beh__r02__0b__neut__nonlure_count',
        'block_0_back_neutface_nonlure_run_2_correct': 'mr_y_tfmri__nback__beh__r02__0b__neut__nonlure__crct_count',
        'block_0_back_neutface_nonlure_run_2_acc': 'mr_y_tfmri__nback__beh__r02__0b__neut__nonlure_acc',
        'block_0_back_neutface_nonlure_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__0b__neut__nonlure_rt',
        'block_0_back_neutface_nonlure_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__0b__neut__nonlure_std',
        'block_0_back_negface_total': 'mr_y_tfmri__nback__beh__0b__neg_count',
        'block_0_back_negface_run_1': 'mr_y_tfmri__nback__beh__r01__0b__neg_count',
        'block_0_back_negface_run_2': 'mr_y_tfmri__nback__beh__r02__0b__neg_count',
        'block_0_back_negface_correct': 'mr_y_tfmri__nback__beh__0b__neg__crct_count',
        'block_0_back_negface_acc': 'mr_y_tfmri__nback__beh__0b__neg_acc',
        'block_0_back_negface_correct_rt': 'mr_y_tfmri__nback__beh__0b__neg_rt',
        'block_0_back_negface_correct_std': 'mr_y_tfmri__nback__beh__0b__neg_std',
        'block_0_back_negface_run_1_correct': 'mr_y_tfmri__nback__beh__r01__0b__neg__crct_count',
        'block_0_back_negface_run_1_acc': 'mr_y_tfmri__nback__beh__r01__0b__neg_acc',
        'block_0_back_negface_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__0b__neg_rt',
        'block_0_back_negface_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__0b__neg_std',
        'block_0_back_negface_run_2_correct': 'mr_y_tfmri__nback__beh__r02__0b__neg__crct_count',
        'block_0_back_negface_run_2_acc': 'mr_y_tfmri__nback__beh__r02__0b__neg_acc',
        'block_0_back_negface_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__0b__neg_rt',
        'block_0_back_negface_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__0b__neg_std',
        'block_0_back_negface_target_total': 'mr_y_tfmri__nback__beh__0b__neg__target_count',
        'block_0_back_negface_target_correct': 'mr_y_tfmri__nback__beh__0b__neg__target__crct_count',
        'block_0_back_negface_target_acc': 'mr_y_tfmri__nback__beh__0b__neg__target_acc',
        'block_0_back_negface_target_rt': 'mr_y_tfmri__nback__beh__0b__neg__target_rt',
        'block_0_back_negface_target_std': 'mr_y_tfmri__nback__beh__0b__neg__target_std',
        'block_0_back_negface_target_run_1_total': 'mr_y_tfmri__nback__beh__r01__0b__neg__target_count',
        'block_0_back_negface_target_run_1_correct': 'mr_y_tfmri__nback__beh__r01__0b__neg__target__crct_count',
        'block_0_back_negface_target_run_1_acc': 'mr_y_tfmri__nback__beh__r01__0b__neg__target_acc',
        'block_0_back_negface_target_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__0b__neg__target_rt',
        'block_0_back_negface_target_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__0b__neg__target_std',
        'block_0_back_negface_target_run_2_total': 'mr_y_tfmri__nback__beh__r02__0b__neg__target_count',
        'block_0_back_negface_target_run_2_correct': 'mr_y_tfmri__nback__beh__r02__0b__neg__target__crct_count',
        'block_0_back_negface_target_run_2_acc': 'mr_y_tfmri__nback__beh__r02__0b__neg__target_acc',
        'block_0_back_negface_target_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__0b__neg__target_rt',
        'block_0_back_negface_target_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__0b__neg__target_std',
        'block_0_back_negface_lure_total': 'mr_y_tfmri__nback__beh__0b__neg__lure_count',
        'block_0_back_negface_lure_correct': 'mr_y_tfmri__nback__beh__0b__neg__lure__crct_count',
        'block_0_back_negface_lure_acc': 'mr_y_tfmri__nback__beh__0b__neg__lure_acc',
        'block_0_back_negface_lure_rt': 'mr_y_tfmri__nback__beh__0b__neg__lure_rt',
        'block_0_back_negface_lure_std': 'mr_y_tfmri__nback__beh__0b__neg__lure_std',
        'block_0_back_negface_lure_run_1_total': 'mr_y_tfmri__nback__beh__r01__0b__neg__lure_count',
        'block_0_back_negface_lure_run_1_correct': 'mr_y_tfmri__nback__beh__r01__0b__neg__lure__crct_count',
        'block_0_back_negface_lure_run_1_acc': 'mr_y_tfmri__nback__beh__r01__0b__neg__lure_acc',
        'block_0_back_negface_lure_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__0b__neg__lure_rt',
        'block_0_back_negface_lure_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__0b__neg__lure_std',
        'block_0_back_negface_lure_run_2_total': 'mr_y_tfmri__nback__beh__r02__0b__neg__lure_count',
        'block_0_back_negface_lure_run_2_correct': 'mr_y_tfmri__nback__beh__r02__0b__neg__lure__crct_count',
        'block_0_back_negface_lure_run_2_acc': 'mr_y_tfmri__nback__beh__r02__0b__neg__lure_acc',
        'block_0_back_negface_lure_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__0b__neg__lure_rt',
        'block_0_back_negface_lure_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__0b__neg__lure_std',
        'block_0_back_negface_nonlure_total': 'mr_y_tfmri__nback__beh__0b__neg__nonlure_count',
        'block_0_back_negface_nonlure_correct': 'mr_y_tfmri__nback__beh__0b__neg__nonlure__crct_count',
        'block_0_back_negface_nonlure_acc': 'mr_y_tfmri__nback__beh__0b__neg__nonlure_acc',
        'block_0_back_negface_nonlure_rt': 'mr_y_tfmri__nback__beh__0b__neg__nonlure_rt',
        'block_0_back_negface_nonlure_std': 'mr_y_tfmri__nback__beh__0b__neg__nonlure_std',
        'block_0_back_negface_nonlure_run_1_total': 'mr_y_tfmri__nback__beh__r01__0b__neg__nonlure_count',
        'block_0_back_negface_nonlure_run_1_correct': 'mr_y_tfmri__nback__beh__r01__0b__neg__nonlure__crct_count',
        'block_0_back_negface_nonlure_run_1_acc': 'mr_y_tfmri__nback__beh__r01__0b__neg__nonlure_acc',
        'block_0_back_negface_nonlure_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__0b__neg__nonlure_rt',
        'block_0_back_negface_nonlure_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__0b__neg__nonlure_std',
        'block_0_back_negface_nonlure_run_2_total': 'mr_y_tfmri__nback__beh__r02__0b__neg__nonlure_count',
        'block_0_back_negface_nonlure_run_2_correct': 'mr_y_tfmri__nback__beh__r02__0b__neg__nonlure__crct_count',
        'block_0_back_negface_nonlure_run_2_acc': 'mr_y_tfmri__nback__beh__r02__0b__neg__nonlure_acc',
        'block_0_back_negface_nonlure_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__0b__neg__nonlure_rt',
        'block_0_back_negface_nonlure_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__0b__neg__nonlure_std',
        'block_0_back_place_total': 'mr_y_tfmri__nback__beh__0b__place_count',
        'block_0_back_place_run_1': 'mr_y_tfmri__nback__beh__r01__0b__place_count',
        'block_0_back_place_run_2': 'mr_y_tfmri__nback__beh__r02__0b__place_count',
        'block_0_back_place_correct': 'mr_y_tfmri__nback__beh__0b__place__crct_count',
        'block_0_back_place_acc': 'mr_y_tfmri__nback__beh__0b__place_acc',
        'block_0_back_place_correct_rt': 'mr_y_tfmri__nback__beh__0b__place_rt',
        'block_0_back_place_correct_std': 'mr_y_tfmri__nback__beh__0b__place_std',
        'block_0_back_place_run_1_correct': 'mr_y_tfmri__nback__beh__r01__0b__place__crct_count',
        'block_0_back_place_run_1_acc': 'mr_y_tfmri__nback__beh__r01__0b__place_acc',
        'block_0_back_place_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__0b__place_rt',
        'block_0_back_place_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__0b__place_std',
        'block_0_back_place_run_2_correct': 'mr_y_tfmri__nback__beh__r02__0b__place__crct_count',
        'block_0_back_place_run_2_acc': 'mr_y_tfmri__nback__beh__r02__0b__place_acc',
        'block_0_back_place_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__0b__place_rt',
        'block_0_back_place_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__0b__place_std',
        'block_0_back_place_target_total': 'mr_y_tfmri__nback__beh__0b__place__target_count',
        'block_0_back_place_target_correct': 'mr_y_tfmri__nback__beh__0b__place__target__crct_count',
        'block_0_back_place_target_acc': 'mr_y_tfmri__nback__beh__0b__place__target_acc',
        'block_0_back_place_target_rt': 'mr_y_tfmri__nback__beh__0b__place__target_rt',
        'block_0_back_place_target_std': 'mr_y_tfmri__nback__beh__0b__place__target_std',
        'block_0_back_place_target_run_1_total': 'mr_y_tfmri__nback__beh__r01__0b__place__target_count',
        'block_0_back_place_target_run_1_correct': 'mr_y_tfmri__nback__beh__r01__0b__place__target__crct_count',
        'block_0_back_place_target_run_1_acc': 'mr_y_tfmri__nback__beh__r01__0b__place__target_acc',
        'block_0_back_place_target_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__0b__place__target_rt',
        'block_0_back_place_target_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__0b__place__target_std',
        'block_0_back_place_target_run_2_total': 'mr_y_tfmri__nback__beh__r02__0b__place__target_count',
        'block_0_back_place_target_run_2_correct': 'mr_y_tfmri__nback__beh__r02__0b__place__target__crct_count',
        'block_0_back_place_target_run_2_acc': 'mr_y_tfmri__nback__beh__r02__0b__place__target_acc',
        'block_0_back_place_target_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__0b__place__target_rt',
        'block_0_back_place_target_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__0b__place__target_std',
        'block_0_back_place_lure_total': 'mr_y_tfmri__nback__beh__0b__place__lure_count',
        'block_0_back_place_lure_correct': 'mr_y_tfmri__nback__beh__0b__place__lure__crct_count',
        'block_0_back_place_lure_acc': 'mr_y_tfmri__nback__beh__0b__place__lure_acc',
        'block_0_back_place_lure_rt': 'mr_y_tfmri__nback__beh__0b__place__lure_rt',
        'block_0_back_place_lure_std': 'mr_y_tfmri__nback__beh__0b__place__lure_std',
        'block_0_back_place_lure_run_1_total': 'mr_y_tfmri__nback__beh__r01__0b__place__lure_count',
        'block_0_back_place_lure_run_1_correct': 'mr_y_tfmri__nback__beh__r01__0b__place__lure__crct_count',
        'block_0_back_place_lure_run_1_acc': 'mr_y_tfmri__nback__beh__r01__0b__place__lure_acc',
        'block_0_back_place_lure_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__0b__place__lure_rt',
        'block_0_back_place_lure_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__0b__place__lure_std',
        'block_0_back_place_lure_run_2_total': 'mr_y_tfmri__nback__beh__r02__0b__place__lure_count',
        'block_0_back_place_lure_run_2_correct': 'mr_y_tfmri__nback__beh__r02__0b__place__lure__crct_count',
        'block_0_back_place_lure_run_2_acc': 'mr_y_tfmri__nback__beh__r02__0b__place__lure_acc',
        'block_0_back_place_lure_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__0b__place__lure_rt',
        'block_0_back_place_lure_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__0b__place__lure_std',
        'block_0_back_place_nonlure_total': 'mr_y_tfmri__nback__beh__0b__place__nonlure_count',
        'block_0_back_place_nonlure_correct': 'mr_y_tfmri__nback__beh__0b__place__nonlure__crct_count',
        'block_0_back_place_nonlure_acc': 'mr_y_tfmri__nback__beh__0b__place__nonlure_acc',
        'block_0_back_place_nonlure_rt': 'mr_y_tfmri__nback__beh__0b__place__nonlure_rt',
        'block_0_back_place_nonlure_std': 'mr_y_tfmri__nback__beh__0b__place__nonlure_std',
        'block_0_back_place_nonlure_run_1_total': 'mr_y_tfmri__nback__beh__r01__0b__place__nonlure_count',
        'block_0_back_place_nonlure_run_1_correct': 'mr_y_tfmri__nback__beh__r01__0b__place__nonlure__crct_count',
        'block_0_back_place_nonlure_run_1_acc': 'mr_y_tfmri__nback__beh__r01__0b__place__nonlure_acc',
        'block_0_back_place_nonlure_run_1_correct_rt': 'mr_y_tfmri__nback__beh__r01__0b__place__nonlure_rt',
        'block_0_back_place_nonlure_run_1_correct_std': 'mr_y_tfmri__nback__beh__r01__0b__place__nonlure_std',
        'block_0_back_place_nonlure_run_2_total': 'mr_y_tfmri__nback__beh__r02__0b__place__nonlure_count',
        'block_0_back_place_nonlure_run_2_correct': 'mr_y_tfmri__nback__beh__r02__0b__place__nonlure__crct_count',
        'block_0_back_place_nonlure_run_2_acc': 'mr_y_tfmri__nback__beh__r02__0b__place__nonlure_acc',
        'block_0_back_place_nonlure_run_2_correct_rt': 'mr_y_tfmri__nback__beh__r02__0b__place__nonlure_rt',
        'block_0_back_place_nonlure_run_2_correct_std': 'mr_y_tfmri__nback__beh__r02__0b__place__nonlure_std',
    },
    'REC': {
        'switch_flag': 'mr_y_tfmri__nback__rec__switch_indicator',
        'newneutface_hr': 'mr_y_tfmri__nback__rec__neut__new_acc',
        'newneutface_fa': 'mr_y_tfmri__nback__rec__neut__new_flsalrm',
        'newplace_hr': 'mr_y_tfmri__nback__rec__place__new_acc',
        'newplace_fa': 'mr_y_tfmri__nback__rec__place__new_flsalrm',
        'newposface_hr': 'mr_y_tfmri__nback__rec__pos__new_acc',
        'newposface_fa': 'mr_y_tfmri__nback__rec__pos__new_flsalrm',
        'newnegface_hr': 'mr_y_tfmri__nback__rec__neg__new_acc',
        'newnegface_fa': 'mr_y_tfmri__nback__rec__neg__new_flsalrm',
        'oldposface_hr': 'mr_y_tfmri__nback__rec__pos_acc',
        'oldposface_fa': 'mr_y_tfmri__nback__rec__pos_flsalrm',
        'oldposface_2_back_hr': 'mr_y_tfmri__nback__rec__2b__pos_acc',
        'oldposface_2_back_fa': 'mr_y_tfmri__nback__rec__2b__pos_flsalrm',
        'oldposface_0_back_hr': 'mr_y_tfmri__nback__rec__0b__pos_acc',
        'oldposface_0_back_fa': 'mr_y_tfmri__nback__rec__0b__pos_flsalrm',
        'oldplace_hr': 'mr_y_tfmri__nback__rec__place_acc',
        'oldplace_fa': 'mr_y_tfmri__nback__rec__place_flsalrm',
        'oldplace_2_back_hr': 'mr_y_tfmri__nback__rec__2b__place_acc',
        'oldplace_2_back_fa': 'mr_y_tfmri__nback__rec__2b__place_flsalrm',
        'oldplace_0_back_hr': 'mr_y_tfmri__nback__rec__0b__place_acc',
        'oldplace_0_back_fa': 'mr_y_tfmri__nback__rec__0b__place_flsalrm',
        'oldnegface_hr': 'mr_y_tfmri__nback__rec__neg_acc',
        'oldnegface_fa': 'mr_y_tfmri__nback__rec__neg_flsalrm',
        'oldnegface_2_back_hr': 'mr_y_tfmri__nback__rec__2b__neg_acc',
        'oldnegface_2_back_fa': 'mr_y_tfmri__nback__rec__2b__neg_flsalrm',
        'oldnegface_0_back_hr': 'mr_y_tfmri__nback__rec__0b__neg_acc',
        'oldnegface_0_back_fa': 'mr_y_tfmri__nback__rec__0b__neg_flsalrm',
        'oldneutface_hr': 'mr_y_tfmri__nback__rec__neut_acc',
        'oldneutface_fa': 'mr_y_tfmri__nback__rec__neut_flsalrm',
        'oldneutface_2_back_hr': 'mr_y_tfmri__nback__rec__2b__neut_acc',
        'oldneutface_2_back_fa': 'mr_y_tfmri__nback__rec__2b__neut_flsalrm',
        'oldneutface_0_back_hr': 'mr_y_tfmri__nback__rec__0b__neut_acc',
        'oldneutface_0_back_fa': 'mr_y_tfmri__nback__rec__0b__neut_flsalrm',
        'posface_pr': 'mr_y_tfmri__nback__rec__pos__corr_acc',
        'posface_br': 'mr_y_tfmri__nback__rec__pos_bias',
        'posface_dprime': 'mr_y_tfmri__nback__rec__pos_dprime',
        'neutface_pr': 'mr_y_tfmri__nback__rec__neut__corr_acc',
        'neutface_br': 'mr_y_tfmri__nback__rec__neut_bias',
        'neutface_dprime': 'mr_y_tfmri__nback__rec__neut_dprime',
        'negface_pr': 'mr_y_tfmri__nback__rec__neg__corr_acc',
        'negface_br': 'mr_y_tfmri__nback__rec__neg_bias',
        'negface_dprime': 'mr_y_tfmri__nback__rec__neg_dprime',
        'place_pr': 'mr_y_tfmri__nback__rec__place__corr_acc',
        'place_br': 'mr_y_tfmri__nback__rec__place_bias',
        'place_dprime': 'mr_y_tfmri__nback__rec__place_dprime',
    },
}



# ── private helpers ────────────────────────────────────────────────────────
def _apply_release_names(behav, task):
    """Return behav with keys renamed to release names; unmapped keys dropped."""
    mapping = _RELEASE_NAMES.get(task.upper(), {})
    return {mapping[k]: v for k, v in behav.items() if k in mapping}



def _to_float(val):
    """Convert TSV cell to float; return NaN for 'n/a' and bad strings."""
    if val is None:
        return np.nan
    s = str(val).strip().lower()
    if s in ('n/a', 'na', 'nan', 'none', ''):
        return np.nan
    try:
        return float(s)
    except ValueError:
        return np.nan


def _to_int(val, default=np.nan):
    f = _to_float(val)
    return int(f) if not math.isnan(f) else default


def _mean(vals):
    """Mean of a numeric iterable, NaN if empty."""
    v = [x for x in vals if not math.isnan(x)]
    return float(np.mean(v)) if v else np.nan


def _sd(vals):
    """Sample standard deviation; 0 for 1 element (MATLAB convention), NaN for 0."""
    v = [x for x in vals if not math.isnan(x)]
    if len(v) >= 2:
        return float(np.std(v, ddof=1))
    return 0.0 if len(v) == 1 else np.nan


def _clip_for_norminv(p):
    """Clip probability to (0.0052, 0.9948) before norminv."""
    lo, hi = 0.0052, 0.9948
    return max(lo, min(hi, p))


def _derive_outstem(fname_json):
    """Derive output stem from JSON filename, stripping _events.json suffix."""
    base = os.path.basename(fname_json)
    base = re.sub(r'_events\.json$', '', base)
    base = re.sub(r'\.json$', '', base)
    return base


def _derive_outstem_tsv(fname_tsv):
    """Derive output stem from TSV filename, stripping the run entity (BIDS
    _run-NN, plus a legacy _scanN, and _desc-invalid) and suffix -- the
    aggregate _behavioral.tsv output has no run entity."""
    base = os.path.basename(fname_tsv)
    base = re.sub(r'_events\.tsv$', '', base)
    base = re.sub(r'\.tsv$', '', base)
    base = re.sub(r'_run-\d+(_desc-invalid)?$', '', base)
    base = re.sub(r'_scan\d+$', '', base)
    return base


def _read_tsv(fname):
    """Read a TSV into a DataFrame; return (df, errcode, errmsg)."""
    try:
        df = pd.read_csv(fname, sep='\t', dtype=str, low_memory=False)
        return df, 0, ''
    except Exception as exc:
        return None, 1, 'cannot read {}: {}'.format(fname, exc)


def _col_float(df, col):
    """Return a float Series for col; NaN for n/a and missing."""
    if col not in df.columns:
        return pd.Series(np.nan, index=df.index)
    return df[col].apply(_to_float)


def _eprime_timing_qc(df):
    """E-Prime metadata and timing QC over the runs used, for MID/SST/NBACK.

    df is the runs used, its run column relabeled to slots 1 and 2. The
    experiment name and versions come from slot 1. The delays are the min
    and max over the runs of each run's constant start_delay,
    start_delay_diff and init_delay_diff, from func/ runs only: a beh/ run
    has no scan, so with no func/ run they are NaN.
    """
    def _first(col, d):
        if col not in d.columns:
            return ''
        v = d[col].dropna().astype(str).str.strip()
        v = v[(v != '') & (v.str.lower() != 'n/a')]
        return v.iloc[0] if len(v) else ''

    d1 = df[df['run'] == '1'] if 'run' in df.columns else df
    expt_name = _first('task_name', d1)
    qc = {
        'eprime_behav_only': int('behavioral' in expt_name.lower()),
        'eprime_expt_name': expt_name,
        'eprime_expt_version': _first('task_version', d1),
        'eprime_runtime_version': _first('runtime_version', d1),
    }
    d_scan = df[~df['_beh'].astype(bool)] if '_beh' in df.columns else df
    runs = sorted(d_scan['run'].unique()) if 'run' in d_scan.columns else []
    for col in ('start_delay', 'start_delay_diff', 'init_delay_diff'):
        vals = [_to_float(_first(col, d_scan[d_scan['run'] == r])) for r in runs]
        vals = [v for v in vals if not math.isnan(v)]
        qc[col + '_min'] = min(vals) if vals else np.nan
        qc[col + '_max'] = max(vals) if vals else np.nan
    return qc


def _write_behav_tsv(behav, fname, verbose):
    """Write a wide-format behavioral dict to a single-row TSV, n/a for
    missing values (NaN, None or empty strings)."""
    behav = dict((k, np.nan if isinstance(v, str) and not v.strip() else v)
                 for k, v in behav.items())
    df = pd.DataFrame([behav], columns=list(behav))
    df.to_csv(fname, sep='\t', index=False, na_rep='n/a')
    if verbose:
        print('wrote {}'.format(fname))


# ── MID behavioral metrics ─────────────────────────────────────────────────

def compute_behav_mid(df):
    """Compute MID behavioral metrics matching abcd_extract_eprime_mid.m.

    Column naming matches MATLAB: {cond}_{phase}_{runlabel}_{metric},
    collapsed groups hiloreward/hiloloss, run labels bothruns/run_1/run_2.
    Returns dict of metric_name → value.
    """
    _CONDITIONS = ['small_reward', 'large_reward', 'small_loss', 'large_loss', 'neutral']
    _COLLAPSED  = [('hiloreward', ['small_reward', 'large_reward']),
                   ('hiloloss',   ['small_loss',   'large_loss'])]
    min_run_resps       = 20
    min_feedback_trials = 4

    df = df[df['trial_type'].str.strip().str.lower() != 'dummy'].copy()
    df['_rt']    = _col_float(df, 'response_time') * 1000
    df['_money'] = _col_float(df, 'money')
    df['_run']   = _col_float(df, 'run').apply(lambda x: int(x) if not math.isnan(x) else 0)
    df['_tt']    = df['trial_type'].str.strip()

    # feedback_window=0 → older MID (no additional RT windows); absent → assume new
    if 'feedback_window' not in df.columns:
        fbflag = True
    else:
        fbflag = bool(_col_float(df, 'feedback_window').eq(1).any())

    # MATLAB pre-initializes both run-1/run-2 slots to NaN and only overwrites
    # the slots for runs that met the trial-count criterion (runs_ok), so a
    # session with only 1 usable run still emits run_2 columns, NaN-filled --
    # rather than omitting them. RUNS is the fixed schema; runs_present is
    # what's actually in the data (Python's own F1 already dropped truncated
    # runs entirely, so runs_present is the MATLAB runs_ok equivalent).
    RUNS = [1, 2]
    runs_present = set(r for r in df['_run'].unique() if r != 0)
    behav = {}

    # ── _pressed/_valid_pressed: response counting, matching
    # check_responses() (response OR premature; response AND NOT
    # premature) -- driven by response_antic, not response_time. Computed
    # early; used by both perform_flag and trialresp counts below.
    early = _col_float(df, 'response_early')
    antic_resp = _col_float(df, 'response_antic')
    df['_pressed'] = (antic_resp == 1) | (early == 1)
    df['_valid_pressed'] = (antic_resp == 1) & (early != 1)

    # ── feedback_flag ──────────────────────────────────────────────────────
    fb_ok = True
    for cond in _CONDITIONS:
        if ((df['_tt'] == cond + '_pos_feedback').sum() < min_feedback_trials or
                (df['_tt'] == cond + '_neg_feedback').sum() < min_feedback_trials):
            fb_ok = False
            break
    behav['feedback_flag'] = 1 if fb_ok else 0

    # ── perform_flag ───────────────────────────────────────────────────────
    # MATLAB uses total trial responses (includes premature) not valid RT
    # count, over BOTH run slots unconditionally -- an absent run contributes
    # 0 responses (empty slice), which correctly fails the threshold, matching
    # MATLAB's check_responses() on event_info already filtered to runs_ok.
    run_total_resps = []
    for r in RUNS:
        run_antic = df[(df['_run'] == r) & df['_tt'].str.endswith('_antic')]
        run_total_resps.append(int(run_antic['_pressed'].sum()))
    behav['perform_flag'] = (
        1 if (run_total_resps and min(run_total_resps) >= min_run_resps) else 0)

    # ── total_numtrials ────────────────────────────────────────────────────
    behav['total_numtrials'] = int(df['_tt'].str.endswith('_antic').sum())
    for r in RUNS:
        behav['total_numtrials_run_{}'.format(r)] = (
            int(((df['_run'] == r) & df['_tt'].str.endswith('_antic')).sum())
            if r in runs_present else np.nan)

    # ── per-condition and collapsed-group metrics ──────────────────────────
    groups = [(c, [c + '_antic'], [c + '_pos_feedback'], [c + '_neg_feedback'])
              for c in _CONDITIONS]
    groups += [(nm, [m + '_antic' for m in mb],
                    [m + '_pos_feedback' for m in mb],
                    [m + '_neg_feedback' for m in mb])
               for nm, mb in _COLLAPSED]

    # MID run labels: bothruns (all runs), then run_1, run_2, ... (fixed RUNS,
    # so an absent run still gets its own key -- an empty slice, matching
    # mean_rt/std_rt/rate's already-NaN-from-empty behavior; numtrials is
    # forced to NaN below since MATLAB pre-initializes it to NaN rather than 0)
    run_variants = [('bothruns', df)] + [('run_{}'.format(r), df[df['_run'] == r])
                                         for r in RUNS]

    for grp_name, antic_tt, pos_tt, neg_tt in groups:
        # anticipation events → MATLAB output phase name "total"
        for rl, rf in run_variants:
            rows    = rf[rf['_tt'].isin(antic_tt)]
            # an RT of exactly 0 is a "not recorded" sentinel, not a real
            # response -- excluded from the mean/std here, but a trial with
            # such a value still counts as a response elsewhere (_pressed
            # above), matching MATLAB's get_rt (RT averaging) vs.
            # check_responses (response counting) split.
            rt_vals = rows.loc[rows['_rt'] != 0, '_rt'].dropna().tolist()
            present = rl == 'bothruns' or int(rl.rsplit('_', 1)[1]) in runs_present
            behav['{}_total_{}_numtrials'.format(grp_name, rl)] = len(rows) if present else np.nan
            behav['{}_total_{}_mean_rt'.format(grp_name, rl)]   = _mean(rt_vals)
            behav['{}_total_{}_std_rt'.format(grp_name, rl)]    = _sd(rt_vals)
        # feedback phases
        for phase, tt, other_tt in [('pos_feedback', pos_tt, neg_tt),
                                     ('neg_feedback', neg_tt, pos_tt)]:
            for rl, rf in run_variants:
                rows  = rf[rf['_tt'].isin(tt)]
                other = rf[rf['_tt'].isin(other_tt)]
                n, n_fb = len(rows), len(rows) + len(other)
                rt_vals = rows.loc[rows['_rt'] != 0, '_rt'].dropna().tolist()
                present = rl == 'bothruns' or int(rl.rsplit('_', 1)[1]) in runs_present
                behav['{}_{}_{}_{}'.format(grp_name, phase, rl, 'numtrials')] = n if present else np.nan
                behav['{}_{}_{}_{}'.format(grp_name, phase, rl, 'rate')] = (
                    n / n_fb if (present and n_fb > 0) else np.nan)
                behav['{}_{}_{}_{}'.format(grp_name, phase, rl, 'mean_rt')] = _mean(rt_vals)
                behav['{}_{}_{}_{}'.format(grp_name, phase, rl, 'std_rt')]  = _sd(rt_vals)

    # ── earnings ───────────────────────────────────────────────────────────
    money_all = df['_money'].dropna()
    behav['earnings_total'] = float(money_all.sum()) if len(money_all) > 0 else np.nan
    for r in RUNS:
        if r not in runs_present:
            behav['earnings_run_{}'.format(r)] = np.nan
            continue
        run_money = df[df['_run'] == r]['_money'].dropna()
        behav['earnings_run_{}'.format(r)] = (
            float(run_money.sum()) if len(run_money) > 0 else np.nan)

    # ── trialresp ──────────────────────────────────────────────────────────
    # Unlike total_numtrials/earnings, MATLAB computes this from
    # check_responses() on event_info already filtered by check_offsets()
    # (drops trials with missing cue onset/offset -- e.g. a trial counted
    # behaviorally despite a truncated raw export, see process_mid's
    # trial_mask), so an absent run naturally yields 0 (not a NaN
    # placeholder) and a timing-truncated trial is excluded here even though
    # it has a real recorded RT and counts toward total_numtrials/accuracy.
    timed = df['_tt'].str.endswith('_antic') & _col_float(df, 'onset').notna()
    for r in RUNS:
        run_antic = df[(df['_run'] == r) & timed]
        behav['trialresp_run_{}_numtrials'.format(r)]       = int(run_antic['_pressed'].sum())
        behav['valid_trialresp_run_{}_numtrials'.format(r)] = int(run_antic['_valid_pressed'].sum())
    all_antic = df[timed]
    behav['trialresp_bothruns_numtrials']      = int(all_antic['_pressed'].sum())
    behav['validtrialresp_bothruns_numtrials'] = int(all_antic['_valid_pressed'].sum())

    # ── old MID format: NaN all anticipation and neg_feedback RT ──────────
    # matches MATLAB fbflag=False logic (abcd_extract_eprime_mid.m lines 763-980)
    if not fbflag:
        for k in list(behav.keys()):
            if ('_mean_rt' in k or '_std_rt' in k) and 'pos_feedback' not in k:
                behav[k] = np.nan

    return behav


# ── NBACK behavioral metrics ───────────────────────────────────────────────

def compute_behav_nback(df):
    """Compute NBACK behavioral metrics matching abcd_extract_eprime_nback.m.

    Uses columns: trial_type, run, response_accuracy, block_type,
                  stim_type, target_type, response_time.
    Returns dict of metric_name → value.
    Column naming matches MATLAB: all_total, block_2_back_total,
    posface_total, block_2_back_posface_total, etc. (no nback_ prefix).
    RT metrics use correct trials only (matching MATLAB correct_rt/correct_std).
    block×stim×target totals use _rt/_std; run-level uses _correct_rt/_correct_std.
    """
    df = df[df['trial_type'].str.strip().str.lower() != 'dummy'].copy()
    df = df[~df['trial_type'].str.strip().str.lower().str.startswith('cue')].copy()

    df['_acc'] = _col_float(df, 'response_accuracy')
    df['_rt']  = _col_float(df, 'response_time') * 1000.0
    df['_bt']  = df.get('block_type',  pd.Series('', index=df.index)).fillna('').str.strip().str.lower()
    df['_st']  = df.get('stim_type',   pd.Series('', index=df.index)).fillna('').str.strip().str.lower()
    df['_tgt'] = df.get('target_type', pd.Series('', index=df.index)).fillna('').str.strip().str.lower()
    df['_run'] = _col_float(df, 'run').apply(lambda x: int(x) if not math.isnan(x) else 0)

    # filter to valid events (block_type, stim_type, target_type all non-empty)
    df = df[(df['_bt'] != '') & (df['_st'] != '') & (df['_tgt'] != '')].copy()

    _TYPES  = ['2_back', '0_back']   # block_type values in TSV (lowercased)
    _BLOCKS = ['2_back', '0_back']   # condnames used in column headers
    _STIMS  = ['posface', 'neutface', 'negface', 'place']
    _TGTS   = ['target', 'lure', 'nonlure']
    RUNS    = [1, 2]

    runs_present = set(int(r) for r in df['_run'].unique() if not math.isnan(float(r)))
    nan = float('nan')
    behav = {}

    def _stats(rows):
        """Return (n_total, n_correct, acc, correct_rt_mean, correct_rt_std)."""
        n = len(rows)
        if n == 0:
            return 0, 0, nan, nan, nan
        cor = rows[rows['_acc'] == 1.0]
        nc  = len(cor)
        acc = nc / n
        rt_vals = cor['_rt'].dropna().tolist()
        return n, nc, acc, _mean(rt_vals), _sd(rt_vals)

    # ── all level (7 metrics) ──────────────────────────────────────────────
    n, nc, ac, rt_m, rt_s = _stats(df)
    behav['all_total'] = n
    for run in RUNS:
        behav['all_run_{}'.format(run)] = (
            len(df[df['_run'] == run]) if run in runs_present else nan)
    behav['all_total_correct']    = nc
    behav['all_total_acc']        = ac
    behav['all_total_correct_rt'] = rt_m
    behav['all_total_correct_std']= rt_s

    # ── block and nested levels ────────────────────────────────────────────
    for bt_lower, bt_cond in zip(_TYPES, _BLOCKS):
        bt_df = df[df['_bt'] == bt_lower]

        # block level
        n, nc, ac, rt_m, rt_s = _stats(bt_df)
        behav['block_{}_total'.format(bt_cond)] = n
        for run in RUNS:
            behav['block_{}_run_{}'.format(bt_cond, run)] = (
                len(bt_df[bt_df['_run'] == run]) if run in runs_present else nan)
        behav['block_{}_total_correct'.format(bt_cond)]    = nc
        behav['block_{}_total_acc'.format(bt_cond)]        = ac
        behav['block_{}_total_correct_rt'.format(bt_cond)] = rt_m
        behav['block_{}_total_correct_std'.format(bt_cond)]= rt_s
        for run in RUNS:
            pfx = 'block_{}_run_{}'.format(bt_cond, run)
            behav['{}_correct'.format(pfx)]    = nan
            behav['{}_acc'.format(pfx)]        = nan
            behav['{}_correct_rt'.format(pfx)] = nan
            behav['{}_correct_std'.format(pfx)]= nan
        for run in sorted(runs_present):
            run_bt = bt_df[bt_df['_run'] == run]
            _, nc_r, ac_r, rt_m_r, rt_s_r = _stats(run_bt)
            pfx = 'block_{}_run_{}'.format(bt_cond, run)
            behav['{}_correct'.format(pfx)]    = nc_r
            behav['{}_acc'.format(pfx)]        = ac_r
            behav['{}_correct_rt'.format(pfx)] = rt_m_r
            behav['{}_correct_std'.format(pfx)]= rt_s_r

        # stim level (inside block loop; uses ALL blocks — last block pass wins,
        # values are identical each pass since ind_stim does not depend on block)
        for st in _STIMS:
            st_df = df[df['_st'] == st]

            n, nc, ac, rt_m, rt_s = _stats(st_df)
            behav['{}_total'.format(st)] = n
            for run in RUNS:
                behav['{}_run_{}'.format(st, run)] = (
                    len(st_df[st_df['_run'] == run]) if run in runs_present else nan)
            behav['{}_correct'.format(st)]    = nc
            behav['{}_acc'.format(st)]        = ac
            behav['{}_correct_rt'.format(st)] = rt_m
            behav['{}_correct_std'.format(st)]= rt_s
            for run in RUNS:
                pfx = '{}_run_{}'.format(st, run)
                behav['{}_correct'.format(pfx)]    = nan
                behav['{}_acc'.format(pfx)]        = nan
                behav['{}_correct_rt'.format(pfx)] = nan
                behav['{}_correct_std'.format(pfx)]= nan
            for run in sorted(runs_present):
                run_st = st_df[st_df['_run'] == run]
                _, nc_r, ac_r, rt_m_r, rt_s_r = _stats(run_st)
                pfx = '{}_run_{}'.format(st, run)
                behav['{}_correct'.format(pfx)]    = nc_r
                behav['{}_acc'.format(pfx)]        = ac_r
                behav['{}_correct_rt'.format(pfx)] = rt_m_r
                behav['{}_correct_std'.format(pfx)]= rt_s_r

            # block×stim level
            bs_df = bt_df[bt_df['_st'] == st]

            n, nc, ac, rt_m, rt_s = _stats(bs_df)
            behav['block_{}_{}_total'.format(bt_cond, st)] = n
            for run in RUNS:
                behav['block_{}_{}_run_{}'.format(bt_cond, st, run)] = (
                    len(bs_df[bs_df['_run'] == run]) if run in runs_present else nan)
            behav['block_{}_{}_correct'.format(bt_cond, st)]    = nc
            behav['block_{}_{}_acc'.format(bt_cond, st)]        = ac
            behav['block_{}_{}_correct_rt'.format(bt_cond, st)] = rt_m
            behav['block_{}_{}_correct_std'.format(bt_cond, st)]= rt_s
            for run in RUNS:
                pfx = 'block_{}_{}_run_{}'.format(bt_cond, st, run)
                behav['{}_correct'.format(pfx)]    = nan
                behav['{}_acc'.format(pfx)]        = nan
                behav['{}_correct_rt'.format(pfx)] = nan
                behav['{}_correct_std'.format(pfx)]= nan
            for run in sorted(runs_present):
                run_bs = bs_df[bs_df['_run'] == run]
                _, nc_r, ac_r, rt_m_r, rt_s_r = _stats(run_bs)
                pfx = 'block_{}_{}_run_{}'.format(bt_cond, st, run)
                behav['{}_correct'.format(pfx)]    = nc_r
                behav['{}_acc'.format(pfx)]        = ac_r
                behav['{}_correct_rt'.format(pfx)] = rt_m_r
                behav['{}_correct_std'.format(pfx)]= rt_s_r

            # block×stim×target level
            for tgt in _TGTS:
                bst_df  = bs_df[bs_df['_tgt'] == tgt]
                bst_pfx = 'block_{}_{}_{}'.format(bt_cond, st, tgt)

                n, nc, ac, rt_m, rt_s = _stats(bst_df)
                # total level: MATLAB uses _rt/_std (not _correct_rt/_correct_std)
                behav['{}_total'.format(bst_pfx)]  = n
                behav['{}_correct'.format(bst_pfx)]= nc
                behav['{}_acc'.format(bst_pfx)]    = ac
                behav['{}_rt'.format(bst_pfx)]     = rt_m
                behav['{}_std'.format(bst_pfx)]    = rt_s
                for run in RUNS:
                    run_pfx = '{}_run_{}'.format(bst_pfx, run)
                    behav['{}_total'.format(run_pfx)]      = nan
                    behav['{}_correct'.format(run_pfx)]    = nan
                    behav['{}_acc'.format(run_pfx)]        = nan
                    behav['{}_correct_rt'.format(run_pfx)] = nan
                    behav['{}_correct_std'.format(run_pfx)]= nan
                for run in sorted(runs_present):
                    run_bst = bst_df[bst_df['_run'] == run]
                    n_r, nc_r, ac_r, rt_m_r, rt_s_r = _stats(run_bst)
                    run_pfx = '{}_run_{}'.format(bst_pfx, run)
                    behav['{}_total'.format(run_pfx)]      = n_r
                    behav['{}_correct'.format(run_pfx)]    = nc_r
                    behav['{}_acc'.format(run_pfx)]        = ac_r
                    behav['{}_correct_rt'.format(run_pfx)] = rt_m_r
                    behav['{}_correct_std'.format(run_pfx)]= rt_s_r

    # ── perform_flag ───────────────────────────────────────────────────────
    acc2 = behav.get('block_2_back_total_acc', nan)
    acc0 = behav.get('block_0_back_total_acc', nan)
    behav['perform_flag'] = (
        1 if (not math.isnan(acc2) and not math.isnan(acc0) and
              acc2 > 0.6 and acc0 > 0.6) else 0)

    return behav


# ── SST behavioral metrics ─────────────────────────────────────────────────

def compute_behav_sst(df):
    """Compute SST behavioral metrics matching abcd_extract_eprime_sst.m.

    Column naming matches MATLAB: CamelCase type names, counts/percent/rt/rt_std
    metrics, total/run_N run labels. Returns dict of metric_name → value.
    """
    _TT_GO   = ['correct_go', 'correctlate_go', 'incorrect_go', 'incorrectlate_go', 'noresp_go']
    _TT_STOP = ['correct_stop', 'incorrect_stop']
    # neither go nor stop: early_stop (STE not recoded as go) and
    #   zero_ssd_response (zero-SSD stop trial with an SSD-window response)
    _TT_OTHER = ['early_stop', 'zero_ssd_response']
    _TT_MAP  = {
        'correct_go':       'CorrectGo',
        'correctlate_go':   'CorrectLateGo',
        'incorrect_go':     'IncorrectGo',
        'incorrectlate_go': 'IncorrectLateGo',
        'noresp_go':        'NoRespGo',
        'correct_stop':     'CorrectStop',
        'incorrect_stop':   'IncorrectStop',
    }
    _NO_RT = {'noresp_go', 'correct_stop'}

    df = df[df['trial_type'].str.strip().str.lower() != 'dummy'].copy()
    df['_tt']    = df['trial_type'].str.strip().str.lower()
    df['_run']   = _col_float(df, 'run').apply(lambda x: int(x) if not math.isnan(x) else 0)
    cgr = (df['response_correct_go'] if 'response_correct_go' in df.columns
           else pd.Series('', index=df.index))
    df['_cgr']   = cgr.fillna('').astype(str).str.strip().str.lower()
    df['_ssd']          = _col_float(df, 'duration_ssd')          * 1000
    df['_ssd_rt']       = _col_float(df, 'response_time_ssd')           * 1000
    df['_go_rt']        = _col_float(df, 'response_time_go_adjusted')   * 1000
    df['_rt']           = _col_float(df, 'response_time')    * 1000
    df['_stop_rt_adj']  = _col_float(df, 'response_time_stop_adjusted') * 1000

    runs = sorted(r for r in df['_run'].unique() if r != 0)
    behav = {}

    # run variants: total first, then run_N (matching MATLAB column order)
    run_variants = [('total', df)] + [('run_{}'.format(r), df[df['_run'] == r]) for r in runs]

    go_rows      = df[df['_tt'].isin(_TT_GO)]
    stop_rows    = df[df['_tt'].isin(_TT_STOP)]
    n_go_total   = len(go_rows)
    n_stop_total = len(stop_rows)

    cgo_rows     = df[df['_tt'] == 'correct_go']
    incgo_rows   = df[df['_tt'] == 'incorrect_go']
    latego_rows  = df[df['_tt'].isin(['correctlate_go', 'incorrectlate_go'])]
    noresp_rows  = df[df['_tt'] == 'noresp_go']
    cs_rows      = df[df['_tt'] == 'correct_stop']
    incstop_rows = df[df['_tt'] == 'incorrect_stop']

    # an RT of exactly 0 is a "not recorded" sentinel (matches MATLAB's
    # mean(all_rt(all_rt~=0)), applied uniformly to go and stop RTs alike)
    cgo_rt_vals     = cgo_rows.loc[cgo_rows['_go_rt'] != 0, '_go_rt'].dropna().tolist()
    incstop_rt_vals = incstop_rows.loc[incstop_rows['_rt'] != 0, '_rt'].dropna().tolist()
    mean_cgo_rt     = _mean(cgo_rt_vals)
    mean_incstop_rt = _mean(incstop_rt_vals)

    # ── perform_flag ───────────────────────────────────────────────────────
    cgo_pct    = len(cgo_rows)    / n_go_total   if n_go_total   > 0 else np.nan
    incgo_pct  = len(incgo_rows)  / n_go_total   if n_go_total   > 0 else np.nan
    latego_pct = len(latego_rows) / n_go_total   if n_go_total   > 0 else np.nan
    noresp_pct = len(noresp_rows) / n_go_total   if n_go_total   > 0 else np.nan
    cs_pct     = len(cs_rows)     / n_stop_total if n_stop_total > 0 else np.nan
    pf = True
    if math.isnan(cgo_pct)    or cgo_pct < 0.6:               pf = False
    if math.isnan(incgo_pct)  or incgo_pct > 0.3:             pf = False
    if math.isnan(latego_pct) or latego_pct > 0.3:            pf = False
    if math.isnan(noresp_pct) or noresp_pct > 0.3:            pf = False
    if math.isnan(cs_pct)     or not (0.2 <= cs_pct <= 0.8):  pf = False
    if n_go_total < 290:                                        pf = False
    if n_stop_total < 40:                                       pf = False
    if (not math.isnan(mean_cgo_rt) and not math.isnan(mean_incstop_rt)
            and mean_cgo_rt < mean_incstop_rt):                pf = False
    behav['perform_flag'] = 1 if pf else 0

    # ── Go_total, Stop_total, Total ────────────────────────────────────────
    behav['Go_total']   = n_go_total
    behav['Stop_total'] = n_stop_total
    # Total also counts the other trials, neither go nor stop (_TT_OTHER);
    # with ste_recode_flag the early_stop trials are go rows
    n_other_total = int(df['_tt'].isin(_TT_OTHER).sum())
    behav['Total']      = n_go_total + n_stop_total + n_other_total

    # ── Go_total_run_N, Stop_total_run_N, Total_run_N (interleaved by run) ─
    for r in runs:
        rf = df[df['_run'] == r]
        n_go   = len(rf[rf['_tt'].isin(_TT_GO)])
        n_stop = len(rf[rf['_tt'].isin(_TT_STOP)])
        behav['Go_total_run_{}'.format(r)]   = n_go
        behav['Stop_total_run_{}'.format(r)] = n_stop
        behav['Total_run_{}'.format(r)] = (
            n_go + n_stop + int(rf['_tt'].isin(_TT_OTHER).sum()))

    # ── per-type: counts, percent, rt (matching MATLAB column order) ───────
    for tt in _TT_GO + _TT_STOP:
        mat    = _TT_MAP[tt]
        is_go  = tt in _TT_GO
        has_rt = tt not in _NO_RT

        # counts (total, run_1, run_2)
        for rl, rf in run_variants:
            behav['{}_counts_{}'.format(mat, rl)] = len(rf[rf['_tt'] == tt])

        # percent (total, run_1, run_2)
        for rl, rf in run_variants:
            n = len(rf[rf['_tt'] == tt])
            denom = (len(rf[rf['_tt'].isin(_TT_GO)])   if is_go
                     else len(rf[rf['_tt'].isin(_TT_STOP)]))
            behav['{}_percent_{}'.format(mat, rl)] = n / denom if denom > 0 else np.nan

        # rt and rt_std interleaved per run variant (total, run_1, run_2)
        if has_rt:
            for rl, rf in run_variants:
                rows = rf[rf['_tt'] == tt]
                rt_col = '_go_rt' if is_go else '_rt'
                rt_vals = rows.loc[rows[rt_col] != 0, rt_col].dropna().tolist()
                behav['{}_rt_{}'.format(mat, rl)]     = _mean(rt_vals)
                behav['{}_rt_std_{}'.format(mat, rl)] = _sd(rt_vals)

    # ── SSD_mean_total: stop trials only; STE excluded regardless of ste_recode_flag ──
    # When ste_recode_flag=False, STE rows have trial_type='early_stop' (not in
    # _TT_STOP), so stop_rows already excludes them.
    # When ste_recode_flag=True, STE rows are correct_go/incorrect_go (not in
    # _TT_STOP), so stop_rows excludes them as well.
    ssd_vals_total = stop_rows['_ssd'].dropna().tolist()
    mean_ssd_total = _mean(ssd_vals_total)
    behav['SSD_mean_total'] = mean_ssd_total

    # ── STE (stop-too-early) trial count/percent ───────────────────────────
    # duration_ssd>0 and response_time_ssd>0, computed directly from the raw per-row columns
    # rather than trial_type=='early_stop', since both are written
    # for STE-recoded go rows too (see _revise_sst_types in
    # cmig_abcd_eprime_to_events.py) -- correct whether or not STE trials
    # were recoded as Go. Percent denominator matches Stop_total (excludes
    # STE itself), mirroring MATLAB's get_behavioral_data().
    for rl, rf in run_variants:
        n_ste  = len(rf[(rf['_ssd'] > 0) & (rf['_ssd_rt'] > 0)])
        n_stop = len(rf[rf['_tt'].isin(_TT_STOP)])
        behav['STE_counts_{}'.format(rl)]  = n_ste
        behav['STE_percent_{}'.format(rl)] = n_ste / n_stop if n_stop > 0 else np.nan

    # ── zero_ssd_response trial count ──────────────────────────────────────
    for rl, rf in run_variants:
        behav['ZeroSSDResp_counts_{}'.format(rl)] = int(
            (rf['_tt'] == 'zero_ssd_response').sum())

    # mean_ssd for SSRT: same as SSD_mean_total (both stop-only now)
    mean_ssd_stop = mean_ssd_total

    # ── SSRT: mean method ─────────────────────────────────────────────────
    ssrt_mean = (mean_cgo_rt - mean_ssd_total
                 if not math.isnan(mean_cgo_rt) and not math.isnan(mean_ssd_total)
                 else np.nan)
    behav['SSRT_mean_total'] = ssrt_mean

    # ── SSRT: integrated method ───────────────────────────────────────────
    # go_rt_max (MATLAB: max([go_trials.go_rt_adjusted]) over all go-type
    # trials before omission fill) must include correctlate_go/
    # incorrectlate_go: their response_time_go_adjusted (go_dur + fix_rt) can exceed
    # regular go RTs, and excluding them understates go_rt_max whenever the
    # SSRT quantile index lands in the omission-filled plateau.
    go_rts_raw  = go_rows['_go_rt'].tolist()
    regular_go  = go_rows[go_rows['_tt'].isin(
        ['correct_go', 'incorrect_go', 'correctlate_go', 'incorrectlate_go'])]
    regular_rts = [_to_float(v) for v in regular_go['_go_rt'].tolist()
                   if not math.isnan(_to_float(v))]
    max_go_rt = max(regular_rts) if regular_rts else np.nan

    go_rts_filled = []
    for v in go_rts_raw:
        f = _to_float(v)
        go_rts_filled.append(max_go_rt if math.isnan(f) else f)
    go_rts_filled = [r for r in go_rts_filled if not math.isnan(r)]

    # n_cs for SSRT uses the correct_stop column (response-field logic),
    # as MATLAB does, not the trial_type count CorrectStop_counts_total uses
    n_correct_stop = int(stop_rows['correct_stop'].apply(
        pd.to_numeric, errors='coerce').fillna(0).sum())
    if n_stop_total > 0 and go_rts_filled:
        prob_stop_fail = 1.0 - n_correct_stop / n_stop_total
        sorted_go = sorted(go_rts_filled)
        n_go_int = len(sorted_go)
        raw_idx  = prob_stop_fail * n_go_int
        if math.ceil(raw_idx) >= n_go_int:
            idxs = [n_go_int - 1]
        else:
            lo = int(math.floor(raw_idx))
            hi = int(math.ceil(raw_idx))
            idxs = list(set([lo, hi]))
        ssrt_int_go = sum(sorted_go[i] for i in idxs) / float(len(idxs))
        # MATLAB SSRT_integrated uses stop-only SSD_mean (mean([stop_trials.ssd_dur]))
        ssrt_int = ssrt_int_go - mean_ssd_stop if not math.isnan(mean_ssd_stop) else np.nan
    else:
        ssrt_int = np.nan
    behav['SSRT_integrated_total'] = ssrt_int

    # ── violatorflag ──────────────────────────────────────────────────────
    # go_rt = mean(response_time_go_adjusted) over responded go trials; stop_rt =
    # mean(response_time_stop_adjusted) over ALL correct_stop==0 trials, not just
    # trial_type=='incorrect_stop', as MATLAB does. Uses its own variable, separate from
    # incstop_rows (IncorrectStop_rt_total, which correctly uses
    # trial_type). violatorflag = stop_rt > go_rt.
    resp_go = go_rows[go_rows['_cgr'].isin(['0', '1'])]
    mean_go_rt_viol   = _mean(resp_go['_go_rt'].dropna().tolist())
    violator_stop_rows = stop_rows[pd.to_numeric(
        stop_rows['correct_stop'], errors='coerce').fillna(1) == 0]
    violator_adj_vals  = violator_stop_rows['_stop_rt_adj'].dropna().tolist()
    mean_violator_adj  = _mean(violator_adj_vals)
    if not math.isnan(mean_violator_adj) and not math.isnan(mean_go_rt_viol):
        behav['violatorflag'] = 1 if mean_violator_adj > mean_go_rt_viol else 0
    elif math.isnan(mean_violator_adj) and not math.isnan(mean_go_rt_viol):
        # No correct_stop==0 trials: MATLAB's NaN > go_rt = False = 0
        behav['violatorflag'] = 0
    else:
        behav['violatorflag'] = np.nan

    # ── glitchflag ────────────────────────────────────────────────────────
    # Stop and zero_ssd_response trials; STE go/early_stop rows excluded regardless
    # of ste_recode_flag. A short response_time_ssd on an STE trial reflects the participant
    # anticipating the stop signal, not a hardware timing glitch.
    glitch_rows = df[df['_tt'].isin(_TT_STOP + ['zero_ssd_response'])]
    ssd_rt_col  = _col_float(glitch_rows, 'response_time_ssd')  * 1000
    ssd_dur_col = _col_float(glitch_rows, 'duration_ssd') * 1000
    glitch = ((ssd_rt_col > 0) & (ssd_rt_col < 50) & (ssd_dur_col <= 50))
    behav['glitchflag'] = 1 if glitch.any() else 0

    # ── zeroSSDcount ──────────────────────────────────────────────────────
    # every trial with duration_ssd 0, zero_ssd_response included, as MATLAB counts
    behav['zeroSSDcount'] = int((df['_ssd'] == 0).sum())

    # ── backfill NaN run columns when fewer than 2 runs present ───────────
    # MATLAB always writes run_1 and run_2 columns (NaN for excluded runs).
    # When only 1 run was valid (e.g., run_2 truncated), add NaN placeholders
    # so output columns match MATLAB schema.
    if len(runs) < 2:
        run2_num = max(runs) + 1 if runs else 2
        for k in list(behav.keys()):
            if '_run_{}_'.format(max(runs) if runs else 1) in k or \
               k.endswith('_run_{}'.format(max(runs) if runs else 1)):
                k2 = k.replace('_run_{}'.format(max(runs) if runs else 1),
                                '_run_{}'.format(run2_num))
                if k2 not in behav:
                    behav[k2] = np.nan

    return behav


# ── REC behavioral metrics ─────────────────────────────────────────────────

def compute_behav_rec(df):
    """Compute REC recognition memory behavioral metrics.

    df contains REC TSV rows (identified by the task tagged from the
    file name). The WM
    condition ('0_back'/'2_back') for each stimulus is read directly from
    the `block_type` column, which Function 1 (process_rec) resolves from
    the raw WM E-Prime file independent of NBACK's own run-validity
    filtering (matching MATLAB's abcd_extract_eprime_rec.m, which re-parses
    the WM file itself rather than depending on NBACK's processed output).

    Uses columns (REC rows): response_accuracy, stim_type, stim, block_type.
    Returns dict of metric_name → value.
    Output columns match MATLAB abcd_extract_eprime_rec.m (excluding metadata).
    """
    # MATLAB recnames order
    _REC_NAMES = ['newneutface', 'newplace', 'newposface', 'newnegface',
                  'oldposface', 'oldplace', 'oldnegface', 'oldneutface']
    _OLD_TYPES = ['oldposface', 'oldplace', 'oldnegface', 'oldneutface']
    # MATLAB stimnames order (for pr/br/dprime)
    _STIM_CATS = ['posface', 'neutface', 'negface', 'place']
    # MATLAB condnames order
    _CONDS = ['2_back', '0_back']

    # split by task (tagged from the file names)
    task_col = df.get('task', pd.Series('', index=df.index)).fillna('').str.upper()
    df_rec = df[task_col == 'REC'].copy()

    df_rec['_acc'] = _col_float(df_rec, 'response_accuracy')
    df_rec['_stype'] = df_rec.get(
        'stim_type', pd.Series('', index=df_rec.index)).fillna('').str.strip().str.lower()
    df_rec['_wm_bt'] = df_rec.get(
        'block_type', pd.Series('n/a', index=df_rec.index)).fillna('n/a').str.strip().str.lower()

    behav = {}

    # ── HR and FA per recognition category (MATLAB recnames loop) ─────────
    hr_by_cat = {}
    fa_by_cat = {}

    for rectype in _REC_NAMES:
        rows = df_rec[df_rec['_stype'] == rectype]
        acc_vals = rows['_acc'].dropna().tolist()
        mean_acc = _mean(acc_vals) if acc_vals else np.nan
        hr = mean_acc
        fa = 1.0 - mean_acc if not math.isnan(mean_acc) else np.nan
        behav['{}_hr'.format(rectype)] = hr
        behav['{}_fa'.format(rectype)] = fa

        if rectype.startswith('old'):
            cat = rectype[3:]
            hr_by_cat[cat] = hr
            # per-WM-condition hr and fa for old stim types
            for cond in _CONDS:
                cond_rows = df_rec[(df_rec['_stype'] == rectype) &
                                   (df_rec['_wm_bt'] == cond)]
                acc_cond = cond_rows['_acc'].dropna().tolist()
                hr_cond = _mean(acc_cond) if acc_cond else np.nan
                fa_cond = 1.0 - hr_cond if not math.isnan(hr_cond) else np.nan
                behav['{}_{}_hr'.format(rectype, cond)] = hr_cond
                behav['{}_{}_fa'.format(rectype, cond)] = fa_cond
        else:
            cat = rectype[3:]
            fa_by_cat[cat] = fa

    # ── Pr, Br, d-prime per stim category (MATLAB stimnames loop) ─────────
    for cat in _STIM_CATS:
        hr = hr_by_cat.get(cat, np.nan)
        fa = fa_by_cat.get(cat, np.nan)
        if math.isnan(hr) or math.isnan(fa):
            continue
        pr = hr - fa
        br = fa / (1.0 - pr) - 0.5 if (1.0 - pr) > 0 else np.nan
        hr_adj = _clip_for_norminv(hr)
        fa_adj = _clip_for_norminv(fa)
        dprime = float(_norm.ppf(hr_adj) - _norm.ppf(fa_adj))
        behav['{}_pr'.format(cat)] = pr
        behav['{}_br'.format(cat)] = br
        behav['{}_dprime'.format(cat)] = dprime

    return behav


# ── private: run-variant helper ────────────────────────────────────────────

def _run_variants(df, run_col, runs):
    """Yield (label, filtered_df) for each run and the 'total' aggregate."""
    for r in runs:
        yield 'run{}'.format(r), df[df[run_col] == r]
    yield 'total', df


# ── public function ────────────────────────────────────────────────────────

def events_to_behavioral(fname_tsv_list, fname_json=None,
                          outdir='.', outstem=None, run_slots=None,
                          forceflag=False, verbose=True):
    """Read events TSV(s) and write behavioral metrics TSV.

    The task is read from each file name (task-<label>), and metadata from
    constant-valued columns in the TSV. fname_json is accepted for backward
    compatibility but is not required.

    Parameters
    ----------
    fname_tsv_list : list of str or str
        One or more _events.tsv paths (one per run), in run order.
    fname_json : str or None
        Path to the _events.json sidecar (optional; used only for the
        output stem).
    outdir : str
        Output directory (created if absent).
    outstem : str or None
        Base output stem; derived from fname_tsv_list[0] if None.
    run_slots : list of int or None
        MID/SST/NBACK: the BIDS run numbers (the run-NN of the file names)
        that fill the release run slots r01 and r02, in that order, e.g.
        [1, 3] when run 2 is invalid, or [2] for a lone run 2 in r01. None
        (default): the runs passed, by ascending run number. At most 2.
    forceflag : bool
        Overwrite existing output files (default False).
    verbose : bool
        Print progress (default True).

    Returns
    -------
    behav_dict : dict
    errcode : int    0=success 1=error
    errmsg : str
    """
    # accept a single filename string as well as a list
    if isinstance(fname_tsv_list, str):
        fname_tsv_list = [fname_tsv_list]

    os.makedirs(outdir, exist_ok=True)

    # ── read and concatenate TSV files ─────────────────────────────────────
    frames = []
    for i, tsv_path in enumerate(fname_tsv_list):
        if not os.path.exists(tsv_path):
            if verbose:
                print('events_to_behavioral: WARNING: not found: {}'.format(tsv_path))
            continue
        df_i, ec, em = _read_tsv(tsv_path)
        if ec:
            if verbose:
                print('events_to_behavioral: WARNING: {}'.format(em))
            continue
        # a beh/ run has no scan, so its timing QC is not reported
        df_i['_beh'] = os.path.basename(os.path.dirname(os.path.abspath(tsv_path))) == 'beh'
        # tag with the BIDS run number from the filename (_run-<NN>_events.tsv),
        # not the file's own run_num column: that column is Function 1's internal
        # per-raw-file index, which resets to 1 for every new raw file, so two
        # different real runs from a multi-file subject-visit could carry the
        # same value and get merged into one run_N bucket below. The filename's
        # BIDS run number is unique per subject-visit by construction.
        m = re.search(r'_run-(\d+)(?:_desc-invalid)?_events\.tsv$',
                      os.path.basename(tsv_path))
        if m:
            df_i['run'] = str(int(m.group(1)))
        elif 'run_num' in df_i.columns:
            df_i['run'] = df_i['run_num']
        else:
            df_i['run'] = str(i + 1)
        # tag with the task from the filename (_task-<label>_), replacing any
        # task column of older files
        m = re.search(r'_task-([A-Za-z0-9]+)_', os.path.basename(tsv_path))
        df_i['task'] = m.group(1).upper() if m else ''
        frames.append(df_i)

    if not frames:
        return {}, 1, 'events_to_behavioral: no TSV files could be read'

    df_all = pd.concat(frames, ignore_index=True, sort=False)

    # exclude dummy, scanner_trigger, and task_init rows
    _special = {'dummy', 'scanner_trigger', 'task_init'}
    df_all = df_all[
        ~df_all['trial_type'].str.strip().str.lower().isin(_special)
    ].copy()

    # ── read task (from the file names) and metadata from TSV columns ──────
    # REC calls pass nBack TSVs + the REC TSV together (for cross-referencing);
    # metadata like switch_flag differs per file, so _meta_col must read from
    # the selected task's own rows, not whichever file happens to appear first
    # in df_all.
    def _meta_col(col, default=''):
        if col not in df_all.columns:
            return default
        df_meta = df_all
        if task and 'task' in df_all.columns:
            df_task = df_all[df_all['task'].str.strip().str.upper() == task]
            if len(df_task) > 0:
                df_meta = df_task
        vals = df_meta[col].dropna()
        vals = vals[vals.str.strip().str.lower() != 'n/a']
        return str(vals.iloc[0]).strip() if len(vals) > 0 else default

    # check all unique task values; prefer REC when present (REC callers pass WM+REC TSVs)
    if 'task' in df_all.columns:
        _task_vals = (df_all['task'].dropna()
                      .str.strip().str.upper()
                      .replace('N/A', float('nan')).dropna().unique().tolist())
        task = 'REC' if 'REC' in _task_vals else (_task_vals[0] if _task_vals else '')
    else:
        task = ''

    if not task:
        return {}, 1, 'events_to_behavioral: no task-<label> in the file names'

    exp_version  = _meta_col('task_version')
    switch_flag  = int(_to_float(_meta_col('switch_flag', '0')) or 0)

    # ── derive outstem ─────────────────────────────────────────────────────
    if outstem is None:
        if fname_json is not None:
            outstem = _derive_outstem(fname_json)
        else:
            outstem = _derive_outstem_tsv(fname_tsv_list[0])

    outstem_path = os.path.join(outdir, outstem)

    if verbose:
        print('events_to_behavioral: task={} runs={} rows={}'.format(
            task, len(frames), len(df_all)))

    # the release schema only ever defines run_1/run_2 columns (MID/SST/
    # NBACK; REC has no run concept), so the runs used are put in slots 1
    # and 2 -- run_slots, else the runs present by run number -- and the run
    # column is relabeled to the slot, which the task metrics read as run_N.
    # More than 2 keeps the first 2, with a warning, rather than inventing
    # non-standard run_3+ release column names. nruns is the runs actually
    # used (<=2).
    n_runs_used = len(frames)
    if task in ('MID', 'SST', 'NBACK') and 'run' in df_all.columns:
        run_vals = pd.to_numeric(df_all['run'], errors='coerce')
        runs_present = sorted(int(r) for r in run_vals.dropna().unique() if r != 0)
        if run_slots is None:
            slot_runs = runs_present
        else:
            slot_runs = [int(r) for r in run_slots]
            missing = [r for r in slot_runs if r not in runs_present]
            if missing:
                print('events_to_behavioral: WARNING: run_slots {} not among '
                      'the runs passed {} for {}'.format(missing, runs_present, outstem))
                slot_runs = [r for r in slot_runs if r in runs_present]
        if len(slot_runs) > 2:
            print('events_to_behavioral: WARNING: {} runs {} for {} -- using '
                  'only the first 2 {}'.format(
                      len(slot_runs), slot_runs, outstem, slot_runs[:2]))
            slot_runs = slot_runs[:2]
        df_all = df_all[run_vals.isin(slot_runs)].copy()
        slot_of = dict((r, i + 1) for i, r in enumerate(slot_runs))
        df_all['run'] = [str(slot_of[int(r)]) for r in
                         pd.to_numeric(df_all['run'], errors='coerce')]
        n_runs_used = len(slot_runs)

    # warn (not an error -- still write the full, correctly-zeroed metrics
    # below) if not a single real event across all runs has a recorded
    # response.
    if 'response_time' in df_all.columns:
        rt_vals = df_all['response_time'].apply(_to_float)
        if not rt_vals.notna().any():
            print('events_to_behavioral: WARNING: no responses recorded in '
                 'any run for {}'.format(outstem))

    # ── dispatch to task-specific metrics ─────────────────────────────────
    if task == 'MID':
        metrics = compute_behav_mid(df_all)
    elif task == 'NBACK':
        metrics = compute_behav_nback(df_all)
    elif task == 'SST':
        metrics = compute_behav_sst(df_all)
    elif task == 'REC':
        metrics = compute_behav_rec(df_all)
    else:
        return {}, 1, 'events_to_behavioral: unsupported task {!r}'.format(task)

    # ── assemble output in MATLAB column order ─────────────────────────────
    # MATLAB schema: SubjID, VisitID, version, switch_flag,
    #   [feedback_flag (MID only),] [perform_flag (MID/SST/NBACK),]
    #   [nruns (MID/SST/NBACK),] ...metric columns...
    behav = {}
    behav['SubjID']      = ''
    behav['VisitID']     = ''
    behav['version']     = exp_version if exp_version else 'UNKNOWN'
    behav['switch_flag'] = switch_flag
    if task == 'MID' and 'feedback_flag' in metrics:
        behav['feedback_flag'] = metrics.pop('feedback_flag')
    if task in ('MID', 'SST', 'NBACK') and 'perform_flag' in metrics:
        behav['perform_flag'] = metrics.pop('perform_flag')
    if task in ('MID', 'SST', 'NBACK'):
        behav['nruns'] = n_runs_used
        behav.update(_eprime_timing_qc(df_all))
    behav.update(metrics)

    # ── apply release names ───────────────────────────────────────────────
    behav = _apply_release_names(behav, task)

    # ── write behavioral TSV ───────────────────────────────────────────────
    out_fname = outstem_path + '_behavioral.tsv'
    if forceflag or not os.path.exists(out_fname):
        _write_behav_tsv(behav, out_fname, verbose=verbose)
    elif verbose:
        print('events_to_behavioral: {} exists, skipping (forceflag=False)'.format(out_fname))

    return behav, 0, ''
