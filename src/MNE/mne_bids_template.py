# Default settings for data processing and analysis.

from collections.abc import Callable, Sequence
from typing import Annotated, Any, Literal

from annotated_types import Ge, Interval, Len, MinLen
from mne import Covariance
from mne_bids import BIDSPath

# from mne_bids_pipeline.typing import (
#     ArbitraryContrast,
#     DigMontageType,
#     FloatArrayLike,
#     PathLike,
# )

# %%
# # General settings - Path
bids_root = "/path/to/BIDS/"

# deriv_root: PathLike | None = None
# """
# The root of the derivatives directory in which the pipeline will store
# the processing results. If `None`, this will be
# `derivatives/mne-bids-pipeline` inside the BIDS root.

# interactive: bool = False
# """
# If True, the steps will provide some interactive elements, such as
# figures. If running the steps from a notebook or Spyder,
# run `%matplotlib qt` in the command line to open the figures in a separate
# window.

# !!! info
#     Enabling interactive mode deactivates parallel processing.
# """

# sessions: list | Literal["all"] = "all"
# """
# The sessions to process. If `'all'`, will process all sessions found in the
# BIDS dataset.
# """
#################################################################################################################
# # General settings - Task/Session/Run
task = "visualoddball"
# task: str = ""
# """
# The task to process.
# """

# runs: Sequence | Literal["all"] = "all"
# """
# The runs to process. If `'all'`, will process all runs found in the
# BIDS dataset.
# """

# Task is rest should be true because we need continuous data
task_is_rest = True
rest_epochs_overlap = 0
rest_epochs_duration = 1
epochs_tmin = 0
baseline = None

# No source estimation needed for us
run_source_estimation = False
#################################################################################################################
# # General settings - Subjects to run
subjects = ["001"]
# subjects: Sequence[str] | Literal["all"] = "all"
# """
# Subjects to analyze. If `'all'`, include all subjects. To only
# include a subset of subjects, pass a list of their identifiers. Even
# if you plan on analyzing only a single subject, pass their identifier
# as a list.

# Please note that if you intend to EXCLUDE only a few subjects, you
# should consider setting `subjects = 'all'` and adding the
# identifiers of the excluded subjects to `exclude_subjects` (see next
# section).

# ???+ example "Example"
#     ```python
#     subjects = 'all'  # Include all subjects.
#     subjects = ['05']  # Only include subject 05.
#     subjects = ['01', '02']  # Only include subjects 01 and 02.
#     ```
# """

# exclude_subjects: Sequence[str] = []
# """
# Specify subjects to exclude from analysis. The MEG empty-room mock-subject
# is automatically excluded from regular analysis.

# ???+ info "Good Practice / Advice"
#     Keep track of the criteria leading you to exclude
#     a participant (e.g. too many movements, missing blocks, aborted experiment,
#     did not understand the instructions, etc, ...)
#     The `emptyroom` subject will be excluded automatically.
# """
#################################################################################################################
# # General setting - types

ch_types = ['eeg']
# ch_types: Annotated[Sequence[Literal["meg", "mag", "grad", "eeg"]], Len(1, 4)] = []
# """
# The channel types to consider.

# ???+ example "Example"
#     ```python
#     # Use EEG channels:
#     ch_types = ['eeg']

#     # Use magnetometer and gradiometer MEG channels:
#     ch_types = ['mag', 'grad']

#     # Use MEG and EEG channels:
#     ch_types = ['meg', 'eeg']
#     ```
# """


data_type = 'eeg'
# data_type: Literal["meg", "eeg"] | None = None
# """
# The BIDS data type.

# For MEG recordings, this will usually be 'meg'; and for EEG, 'eeg'.
# However, if your dataset contains simultaneous recordings of MEG and EEG,
# stored in a single file, you will typically need to set this to 'meg'.
# If `None`, we will assume that the data type matches the channel type.

# ???+ example "Example"
#     The dataset contains simultaneous recordings of MEG and EEG, and we only
#     wish to process the EEG data, which is stored inside the MEG files:

#     ```python
#     ch_types = ['eeg']
#     data_type = 'meg'
#     ```

#     The dataset contains simultaneous recordings of MEG and EEG, and we only
#     wish to process the gradiometer data:

#     ```python
#     ch_types = ['grad']
#     data_type = 'meg'  # or data_type = None
#     ```

#     The dataset contains only EEG data:

#     ```python
#     ch_types = ['eeg']
#     data_type = 'eeg'  # or data_type = None
#     ```
# """
eog_channels = ['HEL', 'HER', 'VER']
# eog_channels: Sequence[str] | None = None
# """
# Specify EOG channels to use, or create virtual EOG channels.
#################################################################################################################

# Set reference
eeg_reference = "average"
# eeg_reference: Literal["average"] | str | Sequence["str"] = "average"
# """
# The EEG reference to use. If `average`, will use the average reference,
# i.e. the average across all channels. If a string, must be the name of a single
# channel. To use multiple channels as reference, set to a list of channel names.

# ???+ example "Example"
#     Use the average reference:
#     ```python
#     eeg_reference = 'average'
#     ```

#     Use the `P9` channel as reference:
#     ```python
#     eeg_reference = 'P9'
#     ```

#     Use the average of the `P9` and `P10` channels as reference:
#     ```python
#     eeg_reference = ['P9', 'P10']
#     ```
# """

# Set montage
eeg_template_montage = "standard_1005"
# eeg_template_montage: str | DigMontageType | None = None
# """
# !!! warning
#     If the data contains channel names that are not part of the template montage, the
#     pipeline run will fail with an error message. You must either pick a different
#     montage or remove those channels via
#     [`drop_channels`][mne_bids_pipeline._config.drop_channels] to continue.


# drop_channels: Sequence[str] = []
# """
# Names of channels to remove from the data. 
# ???+ example "Example"
#     Exclude channels `Fp1` and `Cz` from processing:
#     ```python
#     drop_channels = ['Fp1', 'Cz]
#     ```
# """

analyze_channels = ['eeg']
# analyze_channels: Literal["all", "ch_types"] | Annotated[Sequence["str"], MinLen(1)] = (
#     "ch_types"
# )
# """
# The names of the channels to analyze during ERP/ERF and time-frequency analysis
# steps. 
#################################################################################################################
# """
reader_extra_params = {"units": "uV"}
# reader_extra_params: dict = {}
# """
# Parameters to be passed to `read_raw_bids()` calls when importing raw data.

# ???+ example "Example"
#     Enforce units for EDF files:
#     ```python
#     reader_extra_params = {"units": "uV"}
#     ```
# """

# read_raw_bids_verbose: Literal["error"] | None = None
# """
# Verbosity level to pass to `read_raw_bids(..., verbose=read_raw_bids_verbose)`.
# If you know your dataset will contain files that are not perfectly BIDS
# compliant (e.g., "Did not find any meg.json..."), you can set this to
# `'error'` to suppress warnings emitted by read_raw_bids.
# """
#################################################################################################################

# plot_psd_for_runs: Literal["all"] | Sequence[str] = "all"
# """
# For which runs to add a power spectral density (PSD) plot to the generated
# report. This can take a considerable amount of time if you have many long
# runs. In this case, specify the runs, or pass an empty list to disable raw PSD
# plotting.
# """
random_state = 42
# random_state: int | None = 42
# """
# You can specify the seed of the random number generator (RNG).
# This setting is passed to the ICA algorithm and to the decoding function,
# ensuring reproducible results. Set to `None` to avoid setting the RNG
# to a defined state.
# """
#################################################################################################################
# %%
# # Preprocessing

# ## Break detection

# find_breaks: bool = False
# """
# During an experimental run, the recording might be interrupted by breaks of
# various durations, e.g. to allow the participant to stretch, blink, and swallow
# freely. During these periods, large-scale artifacts are often picked up by the
# recording system. These artifacts can impair certain stages of processing, e.g.
# the peak-detection algorithms we use to find EOG and ECG activity. In some
# cases, even the bad channel detection algorithms might not function optimally.
# It is therefore advisable to mark such break periods for exclusion at early
# processing stages.

# If `True`, try to mark breaks by finding segments of the data where no
# experimental events have occurred. This will then add annotations with the
# description `BAD_break` to the continuous data, causing these segments to be
# ignored in all following processing steps.

# ???+ example "Example"
#     Automatically find break periods, and annotate them as `BAD_break`.
#     ```python
#     find_breaks = True
#     ```

#     Disable break detection.
#     ```python
#     find_breaks = False
#     ```
# """

# min_break_duration: float = 15.0
# """
# The minimal duration (in seconds) of a data segment without any experimental
# events for it to be considered a "break". Note that the minimal duration of the
# generated `BAD_break` annotation will typically be smaller than this, as by
# default, the annotation will not extend across the entire break.
# See [`t_break_annot_start_after_previous_event`][mne_bids_pipeline._config.t_break_annot_start_after_previous_event]
# and [`t_break_annot_stop_before_next_event`][mne_bids_pipeline._config.t_break_annot_stop_before_next_event]
# to control this behavior.

# ???+ example "Example"
#     Periods between two consecutive experimental events must span at least
#     `15` seconds for this period to be considered a "break".
#     ```python
#     min_break_duration = 15.
#     ```
# """  # noqa : E501

# t_break_annot_start_after_previous_event: float = 5.0
# """
# Once a break of at least
# [`min_break_duration`][mne_bids_pipeline._config.min_break_duration]
# seconds has been discovered, we generate a `BAD_break` annotation that does not
# necessarily span the entire break period. Instead, you will typically want to
# start it some time after the last event before the break period, as to not
# unnecessarily discard brain activity immediately following that event.

# This parameter controls how much time (in seconds) should pass after the last
# pre-break event before we start annotating the following segment of the break
# period as bad.

# ???+ example "Example"
#     Once a break period has been detected, add a `BAD_break` annotation to it,
#     starting `5` seconds after the latest pre-break event.
#     ```python
#     t_break_annot_start_after_previous_event = 5.
#     ```

#     Start the `BAD_break` annotation immediately after the last pre-break
#     event.
#     ```python
#     t_break_annot_start_after_previous_event = 0.
#     ```
# """

# t_break_annot_stop_before_next_event: float = 5.0
# """
# Similarly to how
# [`t_break_annot_start_after_previous_event`][mne_bids_pipeline._config.t_break_annot_start_after_previous_event]
# controls the "gap" between beginning of the break period and `BAD_break`
# annotation onset,  this parameter controls how far the annotation should extend
# toward the first experimental event immediately following the break period
# (in seconds). This can help not to waste a post-break trial by marking its
# pre-stimulus period as bad.

# ???+ example "Example"
#     Once a break period has been detected, add a `BAD_break` annotation to it,
#     starting `5` seconds after the latest pre-break event.
#     ```python
#     t_break_annot_start_after_previous_event = 5.
#     ```

#     Start the `BAD_break` annotation immediately after the last pre-break
#     event.
#     ```python
#     t_break_annot_start_after_previous_event = 0.
#     ```
# """
#################################################################################################################
# %%
# ## Bad channel detection
#
# !!! warning
#     This functionality will soon be removed from the pipeline, and
#     will be integrated into MNE-BIDS.
#
# "Bad", i.e. flat and overly noisy channels, can be automatically detected
# using a procedure inspired by the commercial MaxFilter by Elekta. First,
# a copy of the data is low-pass filtered at 40 Hz. Then, channels with
# unusually low variability are flagged as "flat", while channels with
# excessively high variability are flagged as "noisy". Flat and noisy channels
# are marked as "bad" and excluded from subsequent analysis. See
# :func:`mne.preprocssessing.find_bad_channels_maxwell` for more information
# on this procedure. The list of bad channels detected through this procedure
# will be merged with the list of bad channels already present in the dataset,
# if any.

# find_flat_channels_meg: bool = False
# """
# Auto-detect "flat" channels (i.e. those with unusually low variability) and
# mark them as bad.
# """
# find_noisy_channels_meg = True
# find_noisy_channels_meg: bool = False
# """
# Auto-detect "noisy" channels and mark them as bad.
# """
#################################################################################################################

# ## Filtering & resampling

# ### Filtering
#
# It is typically better to set your filtering properties on the raw data so
# as to avoid what we call border (or edge) effects.
#
# If you use this pipeline for evoked responses, you could consider
# a low-pass filter cut-off of h_freq = 40 Hz
# and possibly a high-pass filter cut-off of l_freq = 1 Hz
# so you would preserve only the power in the 1Hz to 40 Hz band.
# Note that highpass filtering is not necessarily recommended as it can
# distort waveforms of evoked components, or simply wash out any low
# frequency that can may contain brain signal. It can also act as
# a replacement for baseline correction in Epochs. See below.
#
# If you use this pipeline for time-frequency analysis, a default filtering
# could be a high-pass filter cut-off of l_freq = 1 Hz
# a low-pass filter cut-off of h_freq = 120 Hz
# so you would preserve only the power in the 1Hz to 120 Hz band.
#
# If you need more fancy analysis, you are already likely past this kind
# of tips! 😇
l_freq = 0.1
# l_freq: float | None = None
# """
# The low-frequency cut-off in the highpass filtering step.
# Keep it `None` if no highpass filtering should be applied.
# """
h_freq = 100.0
# h_freq: float | None = 40.0
# """
# The high-frequency cut-off in the lowpass filtering step.
# Keep it `None` if no lowpass filtering should be applied.
# """

# l_trans_bandwidth: float | Literal["auto"] = "auto"
# """
# Specifies the transition bandwidth of the
# highpass filter. By default it's `'auto'` and uses default MNE
# parameters.
# """

# h_trans_bandwidth: float | Literal["auto"] = "auto"
# """
# Specifies the transition bandwidth of the
# lowpass filter. By default it's `'auto'` and uses default MNE
# parameters.
# """
#################################################################################################################

# # Zapline filter
zapline_fline = 60.0
zapline_iter = True

#################################################################################################################

# ### Resampling
default_resample_sfreq = 250
#
# If you have acquired data with a very high sampling frequency (e.g. 2 kHz)
# you will likely want to downsample to lighten up the size of the files you
# are working with (pragmatics)
# If you are interested in typical analysis (up to 120 Hz) you can typically
# resample your data down to 500 Hz without preventing reliable time-frequency
# exploration of your data.

# raw_resample_sfreq: float | None = None
# """
# Specifies at which sampling frequency the data should be resampled.
# If `None`, then no resampling will be done.

# ???+ example "Example"
#     ```python
#     raw_resample_sfreq = None  # no resampling
#     raw_resample_sfreq = 500  # resample to 500Hz
#     ```
# """
#################################################################################################################

# # ICA
spatial_filter = "ica"
# spatial_filter: Optional[Literal["ssp", "ica"]] = None
# """
# Whether to use a spatial filter to detect and remove artifacts. The BIDS
# Pipeline offers the use of signal-space projection (SSP) and independent
# component analysis (ICA).

# Use `'ssp'` for SSP, `'ica'` for ICA, and `None` if you do not wish to apply
# a spatial filter for artifact removal.

# The Pipeline will try to automatically discover EOG and ECG artifacts. For SSP,
# it will then produce projection vectors that remove ("project out") these
# artifacts from the data. For ICA, the independent components related to
# EOG and ECG activity will be omitted during the signal reconstruction step in
# order to remove the artifacts. The ICA procedure can be configured in various
# ways using the configuration options you can find below.
# """

ica_reject = {'eeg': 800e-6}
# ica_reject: Optional[Union[dict[str, float], Literal["autoreject_local"]]] = None
# """
# Peak-to-peak amplitude limits to exclude epochs from ICA fitting. This allows you to
# remove strong transient artifacts from the epochs used for fitting ICA, which could
# negatively affect ICA performance.

# ???+ info
#     MNE-BIDS-Pipeline will automatically try to detect EOG and ECG artifacts in
#     your data, and remove them. For this to work properly, it is recommended
#     to **not** specify rejection thresholds for EOG and ECG channels here –
#     otherwise, ICA won't be able to "see" these artifacts.

# ???+ info
#     This setting is applied only to the epochs that are used for **fitting** ICA. The
#     goal is to make it easier for ICA to produce a good decomposition. After fitting,
#     ICA is applied to the epochs to be analyzed, usually with one or more components
#     removed (as to remove artifacts). But even after ICA cleaning, some epochs may still
#     contain large-amplitude artifacts. Those epochs can then be rejected by using
#     the [`reject`][mne_bids_pipeline._config.reject] parameter.

# ???+ example "Example"
#     ```python
#     ica_reject = {'grad': 10e-10, 'mag': 20e-12, 'eeg': 400e-6}
#     ica_reject = {'grad': 15e-10}
#     ica_reject = None  # no rejection before fitting ICA
#     ica_reject = "autoreject_global"  # find global (per channel type) PTP thresholds before fitting ICA
#     ica_reject = "autoreject_local"  # find local (per channel) thresholds and repair epochs before fitting ICA
#     ```
# """  # noqa: E501

ica_algorithm = "picard-extended_infomax"
# ica_algorithm: Literal[
#     "picard", "fastica", "extended_infomax", "picard-extended_infomax"
# ] = "picard"
# """
# The ICA algorithm to use. `"picard-extended_infomax"` operates `picard` such that the
# generated ICA decomposition is identical to the one generated by the extended Infomax
# algorithm (but may converge in less time).
# """

ica_l_freq = 1.0
# ica_l_freq: Optional[float] = 1.0
# """
# The cutoff frequency of the high-pass filter to apply before running ICA.
# Using a relatively high cutoff like 1 Hz will remove slow drifts from the
# data, yielding improved ICA results. Must be set to 1 Hz or above.

# Set to `None` to not apply an additional high-pass filter.

# !!! info
#       The filter will be applied to raw data which was already filtered
#       according to the `l_freq` and `h_freq` settings. After filtering, the
#       data will be epoched, and the epochs will be submitted to ICA.

# !!! info
#     The Pipeline will only allow you to perform ICA on data that has been
#     high-pass filtered with a 1 Hz cutoff or higher. This is a conscious,
#     opinionated (but partially data-driven) decision made by the developers.
#     If you have reason to challenge this behavior, please get in touch with
#     us so we can discuss.
# """

ica_max_iterations = 500
# ica_max_iterations: int = 500
# """
# Maximum number of iterations to decompose the data into independent
# components. A low number means to finish earlier, but the consequence is
# that the algorithm may not have finished converging. To ensure
# convergence, pick a high number here (e.g. 3000); yet the algorithm will
# terminate as soon as it determines that is has successfully converged, and
# not necessarily exhaust the maximum number of iterations. Note that the
# default of 200 seems to be sufficient for Picard in many datasets, because
# it converges quicker than the other algorithms; but e.g. for FastICA, this
# limit may be too low to achieve convergence.
# """

ica_n_components = 0.9999999999
# ica_n_components: Optional[Union[float, int]] = None
# """
# MNE conducts ICA as a sort of a two-step procedure: First, a PCA is run
# on the data (trying to exclude zero-valued components in rank-deficient
# data); and in the second step, the principal components are passed
# to the actual ICA. You can select how many of the total principal
# components to pass to ICA – it can be all or just a subset. This determines
# how many independent components to fit, and can be controlled via this
# setting.

# If int, specifies the number of principal components that are passed to the
# ICA algorithm, which will be the number of independent components to
# fit. It must not be greater than the rank of your data (which is typically
# the number of channels, but may be less in some cases).

# If float between 0 and 1, all principal components with cumulative
# explained variance less than the value specified here will be passed to
# ICA.

# If `None` (default), `0.999999` will be used to avoid issues when working with
# rank-deficient data.

# This setting may drastically alter the time required to compute ICA.
# """

# ica_decim: Optional[int] = None
# """
# The decimation parameter to compute ICA. If 5 it means
# that 1 every 5 sample is used by ICA solver. The higher the faster
# it is to run but the less data you have to compute a good ICA. Set to
# `1` or `None` to not perform any decimation.
# """
#################################################################################################################

# # ICLabel
ica_use_icalabel = True
icalabel_include = ["brain", "other"]
ica_use_eog_detection = False
ica_use_ecg_detection = False

#################################################################################################################

# ### Amplitude-based artifact rejection
#
# ???+ info "Good Practice / Advice"
#     Have a look at your raw data and train yourself to detect a blink, a heart
#     beat and an eye movement.
#     You can do a quick average of blink data and check what the amplitude looks
#     like.

reject = {"eeg": 200e-6}
# reject: Optional[
#     Union[dict[str, float], Literal["autoreject_global", "autoreject_local"]]
# ] = None
# """
# Peak-to-peak amplitude limits to mark epochs as bad. This allows you to remove
# epochs with strong transient artifacts.

# !!! info
#       The rejection is performed **after** SSP or ICA, if any of those methods
#       is used. To reject epochs **before** fitting ICA, see the
#       [`ica_reject`][mne_bids_pipeline._config.ica_reject] setting.

# If `None` (default), do not apply artifact rejection.
#################################################################################################################

run_source_estimation = False

#################################################################################################################

mne_log_level = "info"
# log_level: Literal["info", "error"] = "info"
# """
# Set the pipeline logging verbosity.
# """

on_error = "debug"
# on_error: Literal["continue", "abort", "debug"] = "abort"
# """
# Whether to abort processing as soon as an error occurs, continue with all other
# processing steps for as long as possible, or drop you into a debugger in case
# of an error.

# !!! info
#     Enabling debug mode deactivates parallel processing.
# """

