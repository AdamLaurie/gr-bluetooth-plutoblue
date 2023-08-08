/* -*- c++ -*- */
/* 
 * Copyright 2013 Christopher D. Kilgour
 * Copyright 2008, 2009 Dominic Spill, Michael Ossmann
 * Copyright 2007 Dominic Spill
 * Copyright 2005, 2006 Free Software Foundation, Inc.
 *
 * This file is part of gr-bluetooth
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <arrayfire.h>
#include <gnuradio/io_signature.h>
#include "bluetooth/multi_block.h"
#include "bluetooth/packet.h"
#include <gnuradio/filter/firdes.h>
#include <gnuradio/math.h>
#include <stdio.h>
#include <gnuradio/blocks/complex_to_mag_squared.h>

 
namespace gr {
  namespace bluetooth {
    multi_block::multi_block(double sample_rate, double center_freq, double squelch_threshold)
      : gr::sync_block ("bluetooth multi block",
                       gr::io_signature::make (1, 1, sizeof (gr_complex)),
                       gr::io_signature::make (0, 0, 0))
    {
      d_target_snr = squelch_threshold;

      d_cumulative_count = 0;
      d_sample_rate = sample_rate;
      d_center_freq = center_freq;

      /*
       * how many time slots we attempt to decode on each hop:
       * 1 for now, could be as many as 5 plus a little slop
       */
      int slots = 1;
      d_samples_per_symbol = sample_rate / SYMBOL_RATE;
      //FIXME make sure that d_samples_per_symbol >= 2 (requirement of clock_recovery_mm_ff)
      d_samples_per_slot = (int) SYMBOLS_PER_BASIC_RATE_SLOT * d_samples_per_symbol;
      int history_required = (int) slots * d_samples_per_slot;

      /* channel filter coefficients */
      double gain = 1;
      d_channel_filter_width = 500000;
      double transition_width = 300000;
      d_channel_filter = gr::filter::firdes::low_pass( gain,
                                              sample_rate,
                                              d_channel_filter_width,
                                              transition_width,
                                              gr::fft::window::WIN_HANN);
      channel_filter_td=af::array(d_channel_filter.size(),d_channel_filter.data());

      /* noise filter coefficients */
      double n_gain = 1;
      d_noise_filter_width = 22500;
      double n_trans_width = 10000;
      d_noise_filter = gr::filter::firdes::low_pass( n_gain,
                                            sample_rate,
                                            d_noise_filter_width,
                                            n_trans_width,
                                            gr::fft::window::WIN_HANN );
      noise_filter_td=af::array(d_noise_filter.size(),d_noise_filter.data());

      /* we will decimate by the largest integer that results in enough samples per symbol */
      d_ddc_decimation_rate = (int) d_samples_per_symbol / 2;
      double channel_samples_per_symbol = d_samples_per_symbol / d_ddc_decimation_rate;

      set_channels();

      /* fm demodulator */
      d_demod_gain = channel_samples_per_symbol / M_PI_2;

      /* mm_cr variables */
      d_gain_mu = 0.175;
      d_mu = 0.32;
      d_omega_relative_limit = 0.005;
      d_omega = channel_samples_per_symbol;
      d_gain_omega = .25 * d_gain_mu * d_gain_mu;
      d_omega_mid = d_omega;
      d_interp = new gr::filter::mmse_fir_interpolator_ff();
      d_last_sample = 0;
      for (int i=0; i<256; i++) {
	      noise_filtered_samples[i]=0.0;
	      filtered_samples[i]=0.0;
      }
      /* the required history is the slot data + the max of either
         channed DDC + demod, or noise DDC */
      int channel_history = (int) (d_channel_filter.size( ) +
                                   d_ddc_decimation_rate * d_interp->ntaps());
      int noise_history   = (int) d_noise_filter.size( );
      if (channel_history > noise_history) {
        history_required += channel_history;
        d_first_channel_sample = 0;
        d_first_noise_sample   = (channel_history - noise_history);
      }
      else {
        history_required += noise_history;
        d_first_noise_sample   = 0;
        d_first_channel_sample = (noise_history - channel_history);
      }

      printf( "history set to %d samples: channel=%d, noise=%d\n",
              history_required, channel_history, noise_history );

      set_history( history_required );
    }

    static inline float slice(float x)
    {
      return (x < 0) ? -1.0F : 1.0F;
    }

    /* M&M clock recovery, adapted from gr_clock_recovery_mm_ff */
    int
    multi_block::mm_cr(const float *in, int ninput_items, float *out, int noutput_items)
    {
      unsigned int ii = 0; /* input index */
      int          oo = 0; /* output index */
      unsigned int ni = ninput_items - d_interp->ntaps(); /* max input */
      float        mm_val;

      while ((oo < noutput_items) && (ii < ni)) {
        // produce output sample
		//printf("d_mu %3.3f\n", d_mu);
        out[oo]       = d_interp->interpolate( &in[ii], d_mu );
        mm_val        = slice(d_last_sample) * out[oo] - slice(out[oo]) * d_last_sample;
        d_last_sample = out[oo];

        d_omega += d_gain_omega * mm_val;
        d_omega  = d_omega_mid + gr::branchless_clip( d_omega-d_omega_mid,
                                                     d_omega_relative_limit );   // make sure we don't walk away
        d_mu    += d_omega + d_gain_mu * mm_val;

        ii      += (int) floor( d_mu );
        d_mu    -= floor( d_mu );
        oo++;
      }

      /* return number of output items produced */
      return oo;
    }

    /* fm demodulation, taken from gr_quadrature_demod_cf */
    void
    multi_block::demod(const gr_complex *in, float *out, int noutput_items)
    {
      int i;
      gr_complex product;

      for (i = 1; i < noutput_items; i++) {
        gr_complex product = in[i] * conj (in[i-1]);
        out[i] = d_demod_gain * gr::fast_atan2f(imag(product), real(product));
      }
    }

    /* binary slicer, similar to gr_binary_slicer_fb */
    void
    multi_block::slicer(const float *in, char *out, int noutput_items)
    {
      int i;

      for (i = 0; i < noutput_items; i++)
        out[i] = (in[i] < 0) ? 0 : 1;
    }

    int
    multi_block::channel_samples( double                     freq,
                                  gr_vector_const_void_star& in,
                                  gr_vector_void_star&       out,
                                  double&                    energy,
                                  int                        ninput_items )
    {
      if (ninput_items<=d_ddc_decimation_rate) {
	      return 0;
      }
      int ddc_samples = ninput_items - (channel_filter_td.dims(0) - 1) - d_first_channel_sample;
      int ddc_noutput_items=0;
      int classic_chan = abs_freq_channel( freq );
      auto channel_phases=af::constant<float>(0.0f,ddc_samples);
      if (fabs((freq-(d_center_freq))) > 0.00001) {


      	channel_phases=af::array(af::seq((filtered_samples[classic_chan])*((freq)-d_center_freq)*M_PI/(d_sample_rate/2.0),(freq-d_center_freq)*M_PI/(d_sample_rate/2.0)*((ddc_samples+filtered_samples[classic_chan])*2.0),((freq)-d_center_freq)*M_PI/(d_sample_rate/2.0)));

      } 



      channel_phases=channel_phases(af::seq(0.0,ddc_samples-1,1.0));
      auto input_samples_af=af::array(ddc_samples,(af::af_cfloat *)&(((gr_complex*)in[0])[d_first_channel_sample]));
      auto shifted_samples=af::complex(af::cos(channel_phases),af::sin(channel_phases))*input_samples_af;
      auto filtered_samples_real_af=af::fir(channel_filter_td,af::real(shifted_samples));
      auto filtered_samples_imag_af=af::fir(channel_filter_td,af::imag(shifted_samples));
      auto filtered_samples_af=af::complex(filtered_samples_real_af,filtered_samples_imag_af);

      auto decim_mag=af::abs(filtered_samples_af(af::seq(0.0,ddc_samples-1,d_ddc_decimation_rate)));
      auto decim_samples=filtered_samples_af(af::seq(0.0,ddc_samples-1,d_ddc_decimation_rate));
      filtered_samples[classic_chan] += ddc_samples;
      
		// This changes how many iterations it takes to crash... Definitely on to something.
		//printf("ddc_samples: %i\n", ddc_samples);
		//printf("fcs: %i\n", d_first_channel_sample);
		//printf("ddc_noutput_items: %i\n", ddc_noutput_items);
		//gr_vector_void_star ddc_out( 1 );
		//ddc_out[0] = out[0];//malloc(100000);
		//printf("after work %i\n", ddc_noutput_items);
		//printf("past\n");
		//
	ddc_noutput_items=decim_samples.dims(0);
//	printf("noutput_items: %d\n",ddc_noutput_items);
        gr_complex * decim_samples_host=(gr_complex*)decim_samples.host<af::af_cfloat>();
	af::sync();

	float tmpenergy=af::sum<float>(decim_mag);
	tmpenergy=tmpenergy*tmpenergy;
        tmpenergy /= ddc_noutput_items;
	memcpy(out[0],decim_samples_host,ddc_noutput_items*8.0);
	af::freeHost(decim_samples_host);
	energy=(double)(tmpenergy);

		//free(ddc_out[0]);
        //energy /= d_channel_filter_width;

      return ddc_noutput_items;
    }

    int
    multi_block::channel_symbols( gr_vector_const_void_star& in,
                                  char *                     out,
                                  int                        ninput_items )
    {
      /* fm demodulation */
      int demod_noutput_items = ninput_items - 1;
      float demod_out[demod_noutput_items];
      gr_complex *ch_samps = (gr_complex *) in[0];
      demod( ch_samps, demod_out, demod_noutput_items );

      /* clock recovery */
      int cr_ninput_items = demod_noutput_items;
      int noutput_items = cr_ninput_items; // poor estimate but probably safe
      float cr_out[noutput_items];
      noutput_items = mm_cr(demod_out, cr_ninput_items, cr_out, noutput_items);

      /* binary slicer */
      slicer(cr_out, out, noutput_items);

      return noutput_items;
    }

    bool
    multi_block::check_snr( const double               freq,
                            const double               on_channel_energy,
                            double&                    snr,
                            gr_vector_const_void_star& in )
    {
      double off_channel_energy = 0.0;

      int classic_chan = abs_freq_channel( freq );
      auto noise_phases=af::array(af::seq(noise_filtered_samples[classic_chan]*((freq+790000)-d_center_freq)*M_PI/(d_sample_rate/2.0),((freq+790000)-d_center_freq)*M_PI/(d_sample_rate/2.0)*((d_samples_per_slot+noise_filtered_samples[classic_chan])),((freq+790000)-d_center_freq)*M_PI/(d_sample_rate/2.0)));


      int ddc_noutput_items = (int) d_samples_per_slot/d_ddc_decimation_rate;

      noise_phases=noise_phases(af::seq(0.0,d_samples_per_slot-1,1.0));
//      printf("noise_phases: %d\n",noise_phases.dims(0));
      auto input_samples_af=af::array(d_samples_per_slot,(af::af_cfloat *)&(((gr_complex*)in[0])[d_first_noise_sample]));
      auto shifted_samples=af::complex(af::cos(noise_phases),af::sin(noise_phases))*input_samples_af;

      auto filtered_samples_real_af=af::fir(noise_filter_td,af::real(shifted_samples));
      auto filtered_samples_imag_af=af::fir(noise_filter_td,af::imag(shifted_samples));
      auto filtered_samples_af=af::complex(filtered_samples_real_af,filtered_samples_imag_af);
      auto decim_mag=af::abs(filtered_samples_af(af::seq(0.0,d_samples_per_slot-1,d_ddc_decimation_rate)));
	af::sync();
	off_channel_energy=af::sum<float>(decim_mag);
	off_channel_energy=off_channel_energy*off_channel_energy;
        off_channel_energy /= ddc_noutput_items;


        // average mag2 for valley

      snr = 10.0 * log10( on_channel_energy / off_channel_energy );
      if (getenv("SHOW_SNR") && snr >= d_target_snr) {
	printf("%f\n",snr);
      }
      return (snr >= d_target_snr);
    }

    /* add some number of symbols to the block's history requirement */
    void
    multi_block::set_symbol_history(int num_symbols)
    {
      set_history((int) (history() + (num_symbols * d_samples_per_symbol)));
    }

    /* set available channels based on d_center_freq and d_sample_rate */
    void
    multi_block::set_channels()
    {
      /* center frequency described as a fractional channel */
      double center = (d_center_freq - BASE_FREQUENCY) / CHANNEL_WIDTH;
      /* bandwidth in terms of channels */
      double channel_bandwidth = d_sample_rate / CHANNEL_WIDTH;
      /* low edge of our received signal */
      double low_edge = center - (channel_bandwidth / 2);
      /* high edge of our received signal */
      double high_edge = center + (channel_bandwidth / 2);
      /* minimum bandwidth required per channel - ideally 1.0 (1 MHz), but can probably decode with a bit less */
      double min_channel_width = 0.9;

      int low_classic_channel = (int) (low_edge + (min_channel_width / 2) + 1);
      low_classic_channel = (low_classic_channel < 0) ? 0 : low_classic_channel;

      int high_classic_channel = (int) (high_edge - (min_channel_width / 2));
      high_classic_channel = (high_classic_channel > 78) ? 78 : high_classic_channel;

      d_low_freq = channel_abs_freq(low_classic_channel);
      d_high_freq = channel_abs_freq(high_classic_channel);

      for( int ch=low_classic_channel; ch<=high_classic_channel; ch++ ) {
        double freq = channel_abs_freq( ch );

       /* d_channel_ddcs[ch] =
          gr::filter::freq_xlating_fir_filter_ccf::make( d_ddc_decimation_rate,
                                               d_channel_filter,
                                               freq-d_center_freq,
                                               d_sample_rate );
        d_noise_ddcs[ch] =
          gr::filter::freq_xlating_fir_filter_ccf::make( d_ddc_decimation_rate,
                                               d_noise_filter,
                                               freq+790000.0-d_center_freq,
                                               d_sample_rate );*/
      }
    }

    /* returns relative (with respect to d_center_freq) frequency in Hz of given channel */
    double
    multi_block::channel_rel_freq(int channel)
    {
      return channel_abs_freq(channel) - d_center_freq;
    }

    double
    multi_block::channel_abs_freq(int channel)
    {
      return BASE_FREQUENCY + (channel * CHANNEL_WIDTH);
    }

    int
    multi_block::abs_freq_channel(double freq)
    {
      return (int) ((freq - BASE_FREQUENCY) / CHANNEL_WIDTH);
    }

    int multi_block::work (int noutput_items,
                        gr_vector_const_void_star &input_items,
                        gr_vector_void_star &output_items)
    {
        return 0;
    }

  } /* namespace bluetooth */
} /* namespace gr */

