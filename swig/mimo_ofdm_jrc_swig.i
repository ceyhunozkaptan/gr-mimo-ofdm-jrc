/* -*- c++ -*- */

#define MIMO_OFDM_JRC_API

%include "gnuradio.i"           // the common stuff

//load generated python docstrings
%include "mimo_ofdm_jrc_swig_doc.i"

%{
#include "mimo_ofdm_jrc/fft_peak_detect.h"
#include "mimo_ofdm_jrc/frame_detector.h"
#include "mimo_ofdm_jrc/frame_sync.h"
#include "mimo_ofdm_jrc/gui_heatmap_plot.h"
#include "mimo_ofdm_jrc/gui_time_plot.h"
#include "mimo_ofdm_jrc/matrix_transpose.h"
#include "mimo_ofdm_jrc/mimo_ofdm_equalizer.h"
#include "mimo_ofdm_jrc/mimo_ofdm_radar.h"
#include "mimo_ofdm_jrc/mimo_precoder.h"
#include "mimo_ofdm_jrc/moving_avg.h"
#include "mimo_ofdm_jrc/ofdm_cyclic_prefix_remover.h"
#include "mimo_ofdm_jrc/ofdm_frame_generator.h"
#include "mimo_ofdm_jrc/range_angle_estimator.h"
#include "mimo_ofdm_jrc/stream_decoder.h"
#include "mimo_ofdm_jrc/stream_encoder.h"
#include "mimo_ofdm_jrc/target_simulator.h"
#include "mimo_ofdm_jrc/usrp_mimo_trx.h"
#include "mimo_ofdm_jrc/zero_pad.h"
%}

%include "mimo_ofdm_jrc/fft_peak_detect.h"
GR_SWIG_BLOCK_MAGIC2(mimo_ofdm_jrc, fft_peak_detect);
%include "mimo_ofdm_jrc/frame_detector.h"
GR_SWIG_BLOCK_MAGIC2(mimo_ofdm_jrc, frame_detector);
%include "mimo_ofdm_jrc/frame_sync.h"
GR_SWIG_BLOCK_MAGIC2(mimo_ofdm_jrc, frame_sync);
%include "mimo_ofdm_jrc/gui_heatmap_plot.h"
GR_SWIG_BLOCK_MAGIC2(mimo_ofdm_jrc, gui_heatmap_plot);
%include "mimo_ofdm_jrc/gui_time_plot.h"
GR_SWIG_BLOCK_MAGIC2(mimo_ofdm_jrc, gui_time_plot);
%include "mimo_ofdm_jrc/matrix_transpose.h"
GR_SWIG_BLOCK_MAGIC2(mimo_ofdm_jrc, matrix_transpose);
%include "mimo_ofdm_jrc/mimo_ofdm_equalizer.h"
GR_SWIG_BLOCK_MAGIC2(mimo_ofdm_jrc, mimo_ofdm_equalizer);
%include "mimo_ofdm_jrc/mimo_ofdm_radar.h"
GR_SWIG_BLOCK_MAGIC2(mimo_ofdm_jrc, mimo_ofdm_radar);
%include "mimo_ofdm_jrc/mimo_precoder.h"
GR_SWIG_BLOCK_MAGIC2(mimo_ofdm_jrc, mimo_precoder);

%include "mimo_ofdm_jrc/moving_avg.h"
GR_SWIG_BLOCK_MAGIC2(mimo_ofdm_jrc, moving_avg);
%include "mimo_ofdm_jrc/ofdm_cyclic_prefix_remover.h"
GR_SWIG_BLOCK_MAGIC2(mimo_ofdm_jrc, ofdm_cyclic_prefix_remover);
%include "mimo_ofdm_jrc/ofdm_frame_generator.h"
GR_SWIG_BLOCK_MAGIC2(mimo_ofdm_jrc, ofdm_frame_generator);
%include "mimo_ofdm_jrc/range_angle_estimator.h"
GR_SWIG_BLOCK_MAGIC2(mimo_ofdm_jrc, range_angle_estimator);
%include "mimo_ofdm_jrc/stream_decoder.h"
GR_SWIG_BLOCK_MAGIC2(mimo_ofdm_jrc, stream_decoder);
%include "mimo_ofdm_jrc/stream_encoder.h"
GR_SWIG_BLOCK_MAGIC2(mimo_ofdm_jrc, stream_encoder);

%include "mimo_ofdm_jrc/target_simulator.h"
GR_SWIG_BLOCK_MAGIC2(mimo_ofdm_jrc, target_simulator);
%include "mimo_ofdm_jrc/usrp_mimo_trx.h"
GR_SWIG_BLOCK_MAGIC2(mimo_ofdm_jrc, usrp_mimo_trx);
%include "mimo_ofdm_jrc/zero_pad.h"
GR_SWIG_BLOCK_MAGIC2(mimo_ofdm_jrc, zero_pad);

