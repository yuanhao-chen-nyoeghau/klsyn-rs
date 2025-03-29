use crate::*;
use rand::Rng;
use std::f32::consts::PI;

// Filters

pub trait FilterSetPassthrough {
    fn set_passthrough(&mut self);
}
pub trait FilterSetMute {
    fn set_mute(&mut self);
}
pub trait FilterGetTransferFunctionCoefficients {
    /// Returns the polynomial coefficients of the filter transfer function in the z-plane.
    /// The returned array contains the top and bottom coefficients of the rational fraction, ordered in ascending powers.
    fn get_transfer_function_coefficients(&self) -> Fraction;
}
pub trait FilterStep {
    /// Performs a filter step.
    /// # Arguments
    /// * `x` - input signal value
    /// # Returns
    /// * output signal value
    fn step(&mut self, x: f32) -> f32;
}

/// A first-order IIR LP filter.
//
// Formulas:
//  Variables:
//    x = input samples
//    y = output samples
//    a = first filter coefficient
//    b = second filter coefficient, >0 for LP filter, <0 for HP filter
//    f = frequency in Hz
//    w = 2 * PI * f / sampleRate
//    g = gain at frequency f
//  Filter function:
//    y[n] = a * x[n] + b * y[n-1]
//  Transfer function:
//    H(w) = a / ( 1 - b * e^(-jw) )
//  Frequency response:
//    |H(w)| = a / sqrt(1 - 2b * cos(w) + b^2)
//  Gain at DC:
//    |H(0)| = a / sqrt(1 - 2b * cos(0) + b^2)
//           = a / sqrt(1 - 2b + b^2)
//           = a / (1 - b)                                 for b < 1
//  Cutoff frequency for LP filter (frequency with relative gain 0.5, about -3 dB):
//    |H(fCutoff)| = |H(0)| / 2
//    a / sqrt(1 - 2b * cos(w) + b^2) = a / (2 * (1 - b))
//    fCutoff = acos((-3b^2 + 8b - 3) / 2b) * sampleRate / (2 * PI)
//  Determine b for a given gain g at frequency f and |H(0)| = 1:
//    a = 1 - b
//    g = (1 - b) / sqrt(1 - 2b * cos(w) + b^2)
//    g * sqrt(1 - 2b * cos(w) + b^2) = 1 - b
//    g^2 * (1 - 2b * cos(w) + b^2) = 1 - 2b + b^2
//    (g^2 - 1) * b^2  +  2 * (1 - g^2 * cos(w)) * b  +  g^2 - 1  =  0
//    b^2  +  2 * (1 - g^2 * cos(w)) / (g^2 - 1) * b  +  1  =  0
//    Substitute: q = (1 - g^2 * cos(w)) / (1 - g^2)
//    b^2 - 2 * q * b + 1 = 0
//    b = q - sqrt(q^2 - 1)                                or q + sqrt(q^2 - 1)
pub struct LpFilter1 {
    sample_rate: usize,
    /// filter coefficient a
    a: f32,
    /// filter coefficient b
    b: f32,
    /// y[n-1], last output value
    y1: f32,
    passthrough: bool,
    muted: bool,
}
impl LpFilter1 {
    /// # Arguments
    /// * `sample_rate` - sample rate in Hz
    pub fn new(sample_rate: usize) -> Self {
        Self {
            sample_rate,
            a: 0.,
            b: 0.,
            y1: 0.,
            passthrough: true,
            muted: false,
        }
    }
    /// Adjusts the filter parameters without resetting the inner state.
    ///
    /// # Arguments
    /// * `f` - Frequency at which the gain is specified.
    /// * `g` - Gain at frequency f. Between 0 and 1 for LP filter. Greater than 1 for HP filter.
    /// * `extra_gain` - Extra gain factor. This is the resulting DC gain.
    ///   The resulting gain at `f` will be `g * extra_gain`.
    pub fn set(&mut self, f: f32, g: f32, extra_gain: Option<f32>) -> KlResult<()> {
        let extra_gain = extra_gain.unwrap_or(1.);
        if f <= 0.
            || f >= self.sample_rate as f32 / 2.
            || g <= 0.
            || g >= 1.
            || !f.is_finite()
            || !g.is_finite()
            || !extra_gain.is_finite()
        {
            return KlError::err("invalid filter parameters");
        }

        let w = 2. * PI * f / self.sample_rate as f32;
        let gp2 = g.powi(2);
        let q = (1. - gp2 * w.cos()) / (1. - gp2);
        self.b = q - (q.powi(2) - 1.).sqrt();
        self.a = (1. - self.b) * extra_gain;
        self.passthrough = false;
        self.muted = false;

        Ok(())
    }
}
impl FilterSetPassthrough for LpFilter1 {
    fn set_passthrough(&mut self) {
        self.passthrough = true;
        self.muted = false;
        self.y1 = 0.;
    }
}
impl FilterSetMute for LpFilter1 {
    fn set_mute(&mut self) {
        self.passthrough = false;
        self.muted = true;
        self.y1 = 0.;
    }
}
impl FilterGetTransferFunctionCoefficients for LpFilter1 {
    fn get_transfer_function_coefficients(&self) -> Fraction {
        if self.passthrough {
            return fraction!([1], [1]);
        }
        if self.muted {
            return fraction!([0], [1]);
        }
        fraction!([self.a], [1., -self.b])
    }
}
impl FilterStep for LpFilter1 {
    fn step(&mut self, x: f32) -> f32 {
        if self.passthrough {
            return x;
        }
        if self.muted {
            return 0.;
        }
        let y = self.a * x + self.b * self.y1;
        self.y1 = y;
        y
    }
}

/// A Klatt resonator.
/// This is a second order IIR filter.
/// With `f=0` it can also be used as a low-pass filter.
//
// Formulas:
//  Variables:
//    x = input samples
//    y = output samples
//    a/b/c = filter coefficients
//    f = frequency in Hz
//    w = 2 * PI * f / sampleRate
//    f0 = resonator frequency in Hz
//    w0 = 2 * PI * f0 / sampleRate
//    bw = Bandwidth in Hz
//    r = exp(- PI * bw / sampleRate)
//  Filter function:
//    y[n] = a * x[n] + b * y[n-1] + c * y[n-2]
//  Transfer function:
//    H(w) = a / ( 1 - b * e^(-jw) - c * e^(-2jw) )
//  Frequency response:
//    |H(w)| = a / ( sqrt(1 + r^2 - 2 * r * cos(w - w0)) * sqrt(1 + r^2 - 2 * r * cos(w + w0)) )
//  Gain at DC:
//    |H(0)| = a / ( sqrt(1 + r^2 - 2 * r * cos(0 - w0)) * sqrt(1 + r^2 - 2 * r * cos(0 + w0)) )
//           = a / (1 + r^2 - 2 * r * cos(w0))
//           = a / (1 - c - b)
//  Gain at the resonance frequency:
//    |H(f0)| = a / sqrt(1 + r^2 - 2 * r)
//            = a / (1 - r)
#[derive(Clone, Copy)]
pub struct Resonator {
    sample_rate: usize,
    /// filter coefficient a
    a: f32,
    /// filter coefficient b
    b: f32,
    /// filter coefficient c
    c: f32,
    /// y[n-1], last output value
    y1: f32,
    /// y[n-2], second last output value
    y2: f32,
    r: f32,
    passthrough: bool,
    muted: bool,
}
impl Resonator {
    /// # Arguments
    /// * `sample_rate` - sample rate in Hz
    pub fn new(sample_rate: usize) -> Self {
        Self {
            sample_rate,
            a: 0.,
            b: 0.,
            c: 0.,
            y1: 0.,
            y2: 0.,
            r: 0.,
            passthrough: true,
            muted: false,
        }
    }
    /// Adjusts the filter parameters without resetting the inner state.
    ///
    /// # Arguments
    /// * `f` - Frequency of resonator in Hz. May be 0 for LP filtering.
    /// * `bw` - Bandwidth of resonator in Hz.
    /// * `dc_gain` - DC gain level.
    pub fn set(&mut self, f: f32, bw: f32, dc_gain: Option<f32>) -> KlResult<()> {
        let dc_gain = dc_gain.unwrap_or(1.);

        if f < 0.
            || f >= self.sample_rate as f32 / 2.
            || bw <= 0.
            || dc_gain <= 0.
            || !f.is_finite()
            || !bw.is_finite()
            || !dc_gain.is_finite()
        {
            return KlError::err("invalid resonator parameters");
        }

        self.r = (-PI * bw / self.sample_rate as f32).exp();
        let w = 2. * PI * f / self.sample_rate as f32;
        self.c = -self.r.powi(2);
        self.b = 2. * self.r * w.cos();
        self.a = (1. - self.b - self.c) * dc_gain;
        self.passthrough = false;
        self.muted = false;

        Ok(())
    }
    pub fn adjust_impulse_gain(&mut self, new_a: f32) {
        self.a = new_a;
    }
    pub fn adjust_peak_gain(&mut self, peak_gain: f32) -> KlResult<()> {
        if peak_gain <= 0. || !peak_gain.is_finite() {
            return KlError::err("invalid resonator peak gain");
        }
        self.a = peak_gain * (1. - self.r);
        Ok(())
    }
}
impl FilterSetPassthrough for Resonator {
    fn set_passthrough(&mut self) {
        self.passthrough = true;
        self.muted = false;
        self.y1 = 0.;
        self.y2 = 0.;
    }
}
impl FilterSetMute for Resonator {
    fn set_mute(&mut self) {
        self.passthrough = false;
        self.muted = true;
        self.y1 = 0.;
        self.y2 = 0.;
    }
}
impl FilterGetTransferFunctionCoefficients for Resonator {
    fn get_transfer_function_coefficients(&self) -> Fraction {
        if self.passthrough {
            return fraction!([1], [1]);
        }
        if self.muted {
            return fraction!([0], [1]);
        }
        fraction!([self.a], [1., -self.b, -self.c])
    }
}
impl FilterStep for Resonator {
    fn step(&mut self, x: f32) -> f32 {
        if self.passthrough {
            return x;
        }
        if self.muted {
            return 0.;
        }
        let y = self.a * x + self.b * self.y1 + self.c * self.y2;
        self.y2 = self.y1;
        self.y1 = y;
        y
    }
}

/// A Klatt anti-resonator.
/// This is a second order FIR filter.
//
// Formulas:
//  Variables:
//    x = input samples
//    y = output samples
//    a/b/c = filter coefficients
//    f = frequency in Hz
//    w = 2 * PI * f / sampleRate
//  Filter function:
//    y[n] = a * x[n] + b * x[n-1] + c * x[n-2]
//  Transfer function:
//    H(w) = a + b * e^(-jw) + c * e^(-2jw)
struct AntiResonator {
    sample_rate: usize,
    /// filter coefficient a
    a: f32,
    /// filter coefficient b
    b: f32,
    /// filter coefficient c
    c: f32,
    /// x[n-1], last input value
    x1: f32,
    /// x[n-2], second-last input value
    x2: f32,
    passthrough: bool,
    muted: bool,
}
impl AntiResonator {
    /// # Arguments
    /// * `sample_rate` - sample rate in Hz
    pub fn new(sample_rate: usize) -> Self {
        Self {
            sample_rate,
            a: 0.,
            b: 0.,
            c: 0.,
            x1: 0.,
            x2: 0.,
            passthrough: true,
            muted: false,
        }
    }
    /// Adjusts the filter parameters without resetting the inner state.
    ///
    /// # Arguments
    /// * `f` - Frequency of anti-resonator in Hz.
    /// * `bw` - Bandwidth of anti-resonator in Hz.
    pub fn set(&mut self, f: f32, bw: f32) -> KlResult<()> {
        if !valid_freq(f) || f >= self.sample_rate as f32 / 2. || !valid_freq(bw) {
            return KlError::err("invalid anti-resonator parameters");
        }

        let r = (-PI * bw / self.sample_rate as f32).exp();
        let w = 2. * PI * f / self.sample_rate as f32;
        let c0 = -r.powi(2);
        let b0 = 2. * r * w.cos();
        let a0 = 1. - b0 - c0;
        if a0 == 0. {
            self.a = 0.;
            self.b = 0.;
            self.c = 0.;
            return Ok(());
        }

        self.a = 1. / a0;
        self.b = -b0 / a0;
        self.c = -c0 / a0;
        self.passthrough = false;
        self.muted = false;

        Ok(())
    }
}
impl FilterSetPassthrough for AntiResonator {
    fn set_passthrough(&mut self) {
        self.passthrough = true;
        self.muted = false;
        self.x1 = 0.;
        self.x2 = 0.;
    }
}
impl FilterSetMute for AntiResonator {
    fn set_mute(&mut self) {
        self.passthrough = false;
        self.muted = true;
        self.x1 = 0.;
        self.x2 = 0.;
    }
}
impl FilterGetTransferFunctionCoefficients for AntiResonator {
    fn get_transfer_function_coefficients(&self) -> Fraction {
        if self.passthrough {
            return fraction!([1], [1]);
        }
        if self.muted {
            return fraction!([0], [1]);
        }
        fraction!([self.a, self.b, self.c], [1])
    }
}
impl FilterStep for AntiResonator {
    fn step(&mut self, x: f32) -> f32 {
        if self.passthrough {
            return x;
        }
        if self.muted {
            return 0.;
        }
        let y = self.a * x + self.b * self.x1 + self.c * self.x2;
        self.x2 = self.x1;
        self.x1 = x;
        y
    }
}

/// A differencing filter.
/// This is a first-order FIR HP filter.
//
// Problem: The filter curve depends on the sample rate.
// TODO: Compensate the effect of the sample rate.
//
// Formulas:
//  Variables:
//    x = input samples
//    y = output samples
//    f = frequency in Hz
//    w = 2 * PI * f / sampleRate
//  Filter function:
//    y[n] = x[n] - x[n-1]
//  Transfer function:
//    H(w) = 1 - e^(-jw)
//  Frequency response:
//    |H(w)| = sqrt(2 - 2 * cos(w))
struct DifferencingFilter {
    /// x[n-1], last input value
    x1: f32,
}
impl DifferencingFilter {
    pub fn new() -> Self {
        Self { x1: 0. }
    }
}
impl FilterGetTransferFunctionCoefficients for DifferencingFilter {
    fn get_transfer_function_coefficients(&self) -> Fraction {
        fraction!([1., -1.], [1.])
    }
}
impl FilterStep for DifferencingFilter {
    fn step(&mut self, x: f32) -> f32 {
        let y = x - self.x1;
        self.x1 = x;
        y
    }
}

// Sources

trait SourceGetNext {
    fn get_next(&mut self) -> f32;
}

trait SourceStartPeriod {
    /// Starts a new period.
    /// # Arguments
    /// * `open_phase_length` - Duration of the open glottis phase of the F0 period, in samples.
    fn start_period(&mut self, open_phase_length: usize) -> KlResult<()>;
}

// Noise sources

/// Returns a random number within the range -1 .. 1.
fn get_white_noise<R: Rng>(rng: &mut R) -> f32 {
    rng.random_range(-1. ..=1.)
}

/// A low-pass filtered noise source.
pub struct LpNoiseSource<R: Rng> {
    lp_filter: LpFilter1,
    rng: R,
}
impl<R: Rng> LpNoiseSource<R> {
    pub fn new(sample_rate: usize, rng: R) -> KlResult<Self> {
        // The original program logic used a first order LP filter with a filter coefficient
        // of b=0.75 and a sample rate of 10 kHz.
        let old_b = 0.75;
        let old_sample_rate = 10000.;
        // Compute the gain at 1000 Hz with a sample rate of 10 kHz and a DC gain of 1.
        let f = 1000.;
        let g = (1. - old_b)
            / (1. - 2. * old_b * (2. * PI * f / old_sample_rate).cos() + old_b.powi(2)).sqrt();
        // Compensate amplitude for output range -1 .. +1
        let extra_gain = 2.5 * (sample_rate as f32 / 10000.).powf(1. / 3.);
        // Create an LP filter with the same characteristics but with our sampling rate.
        let mut lp_filter = LpFilter1::new(sample_rate);
        lp_filter.set(f, g, Some(extra_gain))?;
        Ok(Self { lp_filter, rng })
    }
}
impl<R: Rng> SourceGetNext for LpNoiseSource<R> {
    fn get_next(&mut self) -> f32 {
        let x = get_white_noise(&mut self.rng);
        self.lp_filter.step(x)
    }
}

// Glottal sources

/// Generates a glottal source signal by LP filtering a pulse train.
pub struct ImpulsiveGlottalSource {
    resonator: Resonator,
    sample_rate: usize,
    position_in_period: usize,
}
impl ImpulsiveGlottalSource {
    pub fn new(sample_rate: usize) -> Self {
        Self {
            resonator: Resonator::new(sample_rate),
            sample_rate,
            position_in_period: 0,
        }
    }
}
impl SourceGetNext for ImpulsiveGlottalSource {
    fn get_next(&mut self) -> f32 {
        if self.resonator.passthrough {
            return 0.;
        }
        let pulse = if self.position_in_period == 1 {
            1.
        } else if self.position_in_period == 2 {
            -1.
        } else {
            0.
        };
        self.position_in_period += 1;
        self.resonator.step(pulse)
    }
}
impl SourceStartPeriod for ImpulsiveGlottalSource {
    fn start_period(&mut self, open_phase_length: usize) -> KlResult<()> {
        if open_phase_length == 0 {
            self.resonator.set_passthrough();
            return Ok(());
        }
        let bw = self.sample_rate as f32 / open_phase_length as f32;
        self.resonator.set(0., bw, None)?;
        self.resonator.adjust_impulse_gain(1.);
        self.position_in_period = 0;
        Ok(())
    }
}

/// Generates a "natural" glottal source signal according to the KLGLOTT88 model.
//
// Formula of the glottal flow: t^2 - t^3
// Formula of the derivative: 2 * t - 3 * t^2
// The derivative is used as the glottal source.
//
// At the end of the open glottal phase there is an abrupt jump from the minimum value to zero.
// This jump is not smoothed in the classic Klatt model. In Praat this "collision phase" is smoothed.
#[derive(Default)]
pub struct NaturalGlottalSource {
    x: f32,
    a: f32,
    b: f32,
    open_phase_length: usize,
    position_in_period: usize,
}
impl NaturalGlottalSource {
    pub fn new() -> Self {
        let mut source = Self::default();
        source.start_period(0).unwrap(); // safe
        source
    }
}
impl SourceGetNext for NaturalGlottalSource {
    fn get_next(&mut self) -> f32 {
        self.position_in_period += 1;
        if self.position_in_period >= self.open_phase_length {
            self.x = 0.;
            return self.x;
        }
        self.a += self.b;
        self.x += self.a;
        self.x
    }
}
impl SourceStartPeriod for NaturalGlottalSource {
    fn start_period(&mut self, open_phase_length: usize) -> KlResult<()> {
        self.open_phase_length = open_phase_length;
        let open_phase_length = open_phase_length as f32;
        self.x = 0.;
        let amplification = 5.;
        self.b = -amplification / open_phase_length.powi(2);
        self.a = -self.b * open_phase_length / 3.;
        self.position_in_period = 0;
        Ok(())
    }
}

pub struct NoiseSource<R: Rng> {
    rng: R,
}
impl<R: Rng> NoiseSource<R> {
    pub fn new(rng: R) -> Self {
        Self { rng }
    }
}
impl<R: Rng> SourceGetNext for NoiseSource<R> {
    fn get_next(&mut self) -> f32 {
        get_white_noise(&mut self.rng)
    }
}
impl<R: Rng> SourceStartPeriod for NoiseSource<R> {
    fn start_period(&mut self, _open_phase_length: usize) -> KlResult<()> {
        Ok(())
    }
}

trait GlottalSource: SourceGetNext + SourceStartPeriod {}
impl GlottalSource for ImpulsiveGlottalSource {}
impl GlottalSource for NaturalGlottalSource {}
impl<R: Rng> GlottalSource for NoiseSource<R> {}

// Utils

/// Modulates the fundamental frequency (F0).
///
/// Sine-wave frequencies of 12.7, 7.1 and 4.7 Hz were chosen so as to ensure
/// a long period before repetition of the perturbation that is introduced.
/// A value of flutter_level = 0.25 results in synthetic vowels with a quite
/// realistic deviation from constant pitch.
///
/// # Arguments
/// * `f0` - Fundamental frequency.
/// * `flutter_level` - Flutter level between 0 and 1.
/// * `time` - Relative signal position in seconds.
/// # Returns
/// * Modulated fundamental frequency.
pub fn perform_frequency_modulation(f0: f32, flutter_level: f32, time: f32) -> f32 {
    if flutter_level <= 0. {
        return f0;
    }
    let w = 2. * PI * time;
    let a = (12.7 * w).sin() + (7.1 * w).sin() + (4.7 * w).sin();
    f0 * (1. + a * flutter_level / 50.)
}

/// Convert a dB value into a linear value.
/// dB values of -99 and below or NaN are converted to 0.
pub fn db_to_lin(db: f32) -> f32 {
    if db <= -99. || db.is_nan() {
        0.
    } else {
        10_f32.powf(db / 20.)
    }
}

fn valid_freq(f: f32) -> bool {
    f > 0. && f.is_finite()
}

// Main logic

pub enum GlottalSourceType {
    Impulsive,
    Natural,
    Noise,
}

pub const MAX_ORAL_FORMANTS: usize = 6;

pub struct MainParms {
    /// Sample rate in Hz.
    pub sample_rate: usize,
    /// Type of glottal source.
    pub glottal_source_type: GlottalSourceType,
}

/// Parameters for a sound frame.
#[derive(PartialEq, PartialOrd)]
pub struct FrameParms {
    /// Frame duration in seconds.
    pub duration: f32,
    /// Fundamental frequency in Hz.
    pub f0: f32,
    /// F0 flutter level, 0 .. 1, typically 0.25.
    pub flutter_level: f32,
    /// Relative length of the open phase of the glottis, 0 .. 1, typically 0.7.
    pub open_phase_ratio: f32,
    /// Breathiness in voicing (turbulence) in dB, positive to amplify or negative to attenuate.
    pub breathiness_db: f32,
    /// Spectral tilt for glottal source in dB. Attenuation at 3 kHz in dB. 0 = no tilt.
    pub tilt_db: f32,
    /// Overall gain (output gain) in dB, positive to amplify, negative to attenuate, `NaN` for automatic gain control (AGC).
    pub gain_db: f32,
    /// RMS level for automatic gain control (AGC), only relevant when `gain_db` is `NaN`.
    pub agc_rms_level: f32,
    /// Nasal formant frequency in Hz, or `NaN`.
    pub nasal_formant_freq: f32,
    /// Nasal formant bandwidth in Hz, or `NaN`.
    pub nasal_formant_bw: f32,
    /// Oral formant frequencies in Hz, or `NaN`.
    pub oral_formant_freq: [f32; MAX_ORAL_FORMANTS],
    /// Oral formant bandwidths in Hz, or `NaN`.
    pub oral_formant_bw: [f32; MAX_ORAL_FORMANTS],

    // Cascade branch:
    /// Whether the cascade branch is enabled.
    pub cascade_enabled: bool,
    /// Voicing amplitude for cascade branch in dB, positive to amplify or negative to attenuate.
    pub cascade_voicing_db: f32,
    /// Aspiration (glottis noise) amplitude for cascade branch in dB, positive to amplify or negative to attenuate.
    pub cascade_aspiration_db: f32,
    /// Amplitude modulation factor for aspiration in cascade branch, 0 = no modulation, 1 = maximum modulation.
    pub cascade_aspiration_mod: f32,
    /// Nasal antiformant frequency in Hz, or `NaN`.
    pub nasal_antiformant_freq: f32,
    /// Nasal antiformant bandwidth in Hz, or `NaN`.
    pub nasal_antiformant_bw: f32,

    // Parallel branch:
    /// Whether the parallel branch is enabled.
    pub parallel_enabled: bool,
    /// Voicing amplitude for parallel branch in dB, positive to amplify or negative to attenuate.
    pub parallel_voicing_db: f32,
    /// Aspiration (glottis noise) amplitude for parallel branch in dB, positive to amplify or negative to attenuate.
    pub parallel_aspiration_db: f32,
    /// Amplitude modulation factor for aspiration in parallel branch, 0 = no modulation, 1 = maximum modulation.
    pub parallel_aspiration_mod: f32,
    /// Frication noise level in dB.
    pub frication_db: f32,
    /// Amplitude modulation factor for frication noise in parallel branch, 0 = no modulation, 1 = maximum modulation.
    pub frication_mod: f32,
    /// Parallel bypass level in dB, used to bypass differentiated glottal and frication signals around resonators F2 to F6.
    pub parallel_bypass_db: f32,
    /// Nasal formant level in dB.
    pub nasal_formant_db: f32,
    /// Oral formant levels in dB, or `NaN`.
    pub oral_formant_db: [f32; MAX_ORAL_FORMANTS],
}

/// Variables of the currently active frame.
#[derive(Default)]
struct FrameState {
    /// Linear breathiness level.
    breathiness_lin: f32,
    /// Linear overall gain.
    gain_lin: f32,

    // Cascade branch:
    /// Linear voicing amplitude for cascade branch.
    cascade_voicing_lin: f32,
    /// Linear aspiration amplitude for cascade branch.
    cascade_aspiration_lin: f32,

    // Parallel branch:
    /// Linear voicing amplitude for parallel branch.
    parallel_voicing_lin: f32,
    /// Linear aspiration amplitude for parallel branch.
    parallel_aspiration_lin: f32,
    /// Linear frication noise level.
    frication_lin: f32,
    /// Linear parallel bypass level.
    parallel_bypass_lin: f32,
}

/// Variables of the currently active F0 period (aka glottal period).
#[derive(Default)]
struct PeriodState {
    /// Modulated fundamental frequency for this period, in Hz, or 0.
    f0: f32,
    /// Period length in samples.
    period_length: usize,
    /// Open glottis phase length in samples.
    open_phase_length: usize,
    /// Current sample position within F0 period.
    position_in_period: usize,
    /// LP filtered noise.
    #[allow(dead_code)]
    lp_noise: f32,
}

/// Sound generator controller.
pub struct Generator<'a, R: Rng + Clone + 'static> {
    /// Main parameters.
    m_parms: &'a MainParms,
    /// Currently active frame parameters.
    f_parms: Option<&'a FrameParms>,
    /// New frame parameters for start of next F0 period.
    new_f_parms: Option<&'a FrameParms>,
    /// Frame variables.
    f_state: FrameState,
    /// F0 period state variables.
    p_state: Option<PeriodState>,
    /// Current absolute sample position.
    abs_position: usize,
    /// Spectral tilt filter.
    tilt_filter: LpFilter1,
    /// Output low-pass filter.
    output_lp_filter: Resonator,
    /// Random value for flutter time offset.
    flutter_time_offset: usize,

    // Glottal sources:
    glottal_source: Box<dyn GlottalSource>,

    // Noise sources:
    // (We use independent noise sources to avoid cancellation effects of correlated signals.)
    /// Noise source for aspiration in cascade branch.
    aspiration_source_casc: LpNoiseSource<R>,
    /// Noise source for aspiration in parallel branch.
    aspiration_source_par: LpNoiseSource<R>,
    /// Noise source for frication in parallel branch.
    frication_source_par: LpNoiseSource<R>,

    // Cascade branch variables:
    /// Nasal formant filter for cascade branch.
    nasal_formant_casc: Resonator,
    /// Nasal antiformant filter for cascade branch.
    nasal_antiformant_casc: AntiResonator,
    /// Oral formant filters for cascade branch.
    oral_formant_casc: [Resonator; MAX_ORAL_FORMANTS],

    // Parallel branch variables:
    /// Nasal formant filter for parallel branch.
    nasal_formant_par: Resonator,
    /// Oral formant filters for parallel branch.
    oral_formant_par: [Resonator; MAX_ORAL_FORMANTS],
    /// Differencing filter for the parallel branch.
    differencing_filter_par: DifferencingFilter,

    rng: R,
}
impl<'a, R: Rng + Clone + 'static> Generator<'a, R> {
    /// Creates a new generator.
    /// # Arguments
    /// * `m_parms` - Main parameters.
    /// * `rng` - Random number generator.
    /// # Returns
    /// * `Ok(Generator)` on success.
    /// * `Err(KlError)` on error.
    pub fn new(m_parms: &'a MainParms, rng: R) -> KlResult<Self> {
        let mut generator = Self {
            m_parms,
            f_parms: None,
            new_f_parms: None,
            f_state: FrameState::default(),
            p_state: None,
            abs_position: 0,
            tilt_filter: LpFilter1::new(m_parms.sample_rate),
            flutter_time_offset: rng.clone().random_range(0..1000),
            glottal_source: match m_parms.glottal_source_type {
                GlottalSourceType::Impulsive => {
                    Box::new(ImpulsiveGlottalSource::new(m_parms.sample_rate))
                }
                GlottalSourceType::Natural => Box::new(NaturalGlottalSource::new()),
                GlottalSourceType::Noise => Box::new(NoiseSource::new(rng.clone())),
            },
            output_lp_filter: Resonator::new(m_parms.sample_rate),
            aspiration_source_casc: LpNoiseSource::new(m_parms.sample_rate, rng.clone())?,
            aspiration_source_par: LpNoiseSource::new(m_parms.sample_rate, rng.clone())?,
            frication_source_par: LpNoiseSource::new(m_parms.sample_rate, rng.clone())?,
            nasal_formant_casc: Resonator::new(m_parms.sample_rate),
            nasal_antiformant_casc: AntiResonator::new(m_parms.sample_rate),
            oral_formant_casc: [Resonator::new(m_parms.sample_rate); MAX_ORAL_FORMANTS],
            nasal_formant_par: Resonator::new(m_parms.sample_rate),
            oral_formant_par: [Resonator::new(m_parms.sample_rate); MAX_ORAL_FORMANTS],
            differencing_filter_par: DifferencingFilter::new(),
            rng,
        };

        generator
            .output_lp_filter
            .set(0., m_parms.sample_rate as f32 / 2., None)?;

        Ok(generator)
    }

    /// Generates a frame of the sound.
    /// The length of the frame is specified by `out_buf.len()` and `f_parms.duration` is ignored.
    /// # Arguments
    /// * `f_parms` - Frame parameters.
    /// * `out_buf` - Output buffer.
    /// # Returns
    /// * `Ok(())` on success.
    /// * `Err(KlError)` on error.
    pub fn generate_frame(&mut self, f_parms: &'a FrameParms, out_buf: &mut [f32]) -> KlResult<()> {
        if let Some(self_f_parms) = self.f_parms {
            if f_parms == self_f_parms {
                return KlError::err("FrameParms structure must not be re-used.");
            }
        }
        self.new_f_parms = Some(f_parms);
        for out_pos in 0..out_buf.len() {
            if self.p_state.is_none()
                || self.p_state.as_ref().unwrap().position_in_period
                    >= self.p_state.as_ref().unwrap().period_length
            {
                self.start_new_period()?;
            }
            out_buf[out_pos] = self.compute_next_output_signal_sample();
            let p_state = self.p_state.as_mut().unwrap(); // safe
            p_state.position_in_period += 1;
            self.abs_position += 1;
        }
        if f_parms.gain_db.is_nan() {
            // automatic gain control (AGC)
            adjust_signal_gain(out_buf, f_parms.agc_rms_level);
        }
        Ok(())
    }

    fn compute_next_output_signal_sample(&mut self) -> f32 {
        let f_parms = self.f_parms.unwrap();
        let p_state = self.p_state.as_ref().unwrap();
        let mut voice = self.glottal_source.get_next();
        voice = self.tilt_filter.step(voice); // apply spectral tilt
        if p_state.position_in_period < p_state.open_phase_length {
            // if within glottal open phase
            voice += get_white_noise(&mut self.rng) * self.f_state.breathiness_lin;
            // add breathiness (turbulence)
        }
        let cascade_out = if f_parms.cascade_enabled {
            self.compute_cascade_branch(voice)
        } else {
            0.
        };
        let parallel_out = if f_parms.parallel_enabled {
            self.compute_parallel_branch(voice)
        } else {
            0.
        };
        let mut out = cascade_out + parallel_out;
        out = self.output_lp_filter.step(out);
        out *= self.f_state.gain_lin;
        out
    }

    fn compute_cascade_branch(&mut self, voice: f32) -> f32 {
        let f_parms = self.f_parms.unwrap();
        let p_state = self.p_state.as_ref().unwrap();
        let cascade_voice = voice * self.f_state.cascade_voicing_lin;
        let current_aspiration_mod = if p_state.position_in_period >= p_state.period_length / 2 {
            f_parms.cascade_aspiration_mod
        } else {
            0.
        };
        let aspiration = self.aspiration_source_casc.get_next()
            * self.f_state.cascade_aspiration_lin
            * (1. - current_aspiration_mod);
        let mut v = cascade_voice + aspiration;
        v = self.nasal_antiformant_casc.step(v);
        v = self.nasal_formant_casc.step(v);
        for i in 0..MAX_ORAL_FORMANTS {
            v = self.oral_formant_casc[i].step(v);
        }
        v
    }

    fn compute_parallel_branch(&mut self, voice: f32) -> f32 {
        let f_parms = self.f_parms.unwrap();
        let p_state = self.p_state.as_ref().unwrap();
        let parallel_voice = voice * self.f_state.parallel_voicing_lin;
        let current_aspiration_mod = if p_state.position_in_period >= p_state.period_length / 2 {
            f_parms.parallel_aspiration_mod
        } else {
            0.
        };
        let aspiration = self.aspiration_source_par.get_next()
            * self.f_state.parallel_aspiration_lin
            * (1. - current_aspiration_mod);
        let source = parallel_voice + aspiration;
        let source_difference = self.differencing_filter_par.step(source);
        // Klatt (1980) states: "... using a first difference calculation to remove low-frequency energy from
        // the higher formants; this energy would otherwise distort the spectrum in the region of F1 during
        // the synthesis of some vowels."
        // A differencing filter is applied for H2 to H6 and the bypass.
        // A better solution would probably be to use real band-pass filters instead of resonators for the formants
        // in the parallel branch. Then this differencing filter would not be necessary to protect the low frequencies
        // of the low formants.
        let current_frication_mod = if p_state.position_in_period >= p_state.period_length / 2 {
            f_parms.frication_mod
        } else {
            0.
        };
        let frication_noise = self.frication_source_par.get_next()
            * self.f_state.frication_lin
            * (1. - current_frication_mod);
        let source2 = source_difference + frication_noise;
        let mut v = 0.;
        v += self.nasal_formant_par.step(source); // nasal formant is directly applied to source
        v += self.oral_formant_par[0].step(source); // F1 is directly applied to source
        for i in 1..MAX_ORAL_FORMANTS {
            let alternating_sign = if i % 2 == 0 { 1. } else { -1. };
            v += alternating_sign * self.oral_formant_par[i].step(source2);
        }
        v += self.f_state.parallel_bypass_lin * source2;
        v
    }

    /// Starts a new F0 period.
    fn start_new_period(&mut self) -> KlResult<()> {
        if let Some(f_parms) = self.new_f_parms {
            self.f_parms = Some(f_parms);
            self.new_f_parms = None;
            self.start_using_new_frame_parameters()?;
        }
        if self.p_state.is_none() {
            self.p_state = Some(PeriodState::default());
        }
        let p_state = self.p_state.as_mut().unwrap(); // safe
        let m_parms = self.m_parms;
        let f_parms = self.f_parms.unwrap();
        let flutter_time =
            self.abs_position as f32 / m_parms.sample_rate as f32 + self.flutter_time_offset as f32;
        p_state.f0 = perform_frequency_modulation(f_parms.f0, f_parms.flutter_level, flutter_time);
        p_state.period_length = if p_state.f0 > 0. {
            (m_parms.sample_rate as f32 / p_state.f0).round() as usize
        } else {
            1
        };
        p_state.open_phase_length = if p_state.period_length > 1 {
            (p_state.period_length as f32 * f_parms.open_phase_ratio).round() as usize
        } else {
            0
        };
        p_state.position_in_period = 0;
        self.start_glottal_source_period()?;
        Ok(())
    }

    fn start_using_new_frame_parameters(&mut self) -> KlResult<()> {
        let m_parms = self.m_parms;
        let f_parms = self.f_parms.unwrap();
        self.f_state.breathiness_lin = db_to_lin(f_parms.breathiness_db);
        self.f_state.gain_lin = db_to_lin(if f_parms.gain_db.is_nan() {
            0.
        } else {
            f_parms.gain_db
        });
        set_tilt_filter(&mut self.tilt_filter, f_parms.tilt_db)?;

        // Adjust cascade branch:
        self.f_state.cascade_voicing_lin = db_to_lin(f_parms.cascade_voicing_db);
        self.f_state.cascade_aspiration_lin = db_to_lin(f_parms.cascade_aspiration_db);
        set_nasal_formant_casc(&mut self.nasal_formant_casc, f_parms)?;
        set_nasal_antiformant_casc(&mut self.nasal_antiformant_casc, f_parms)?;
        for i in 0..MAX_ORAL_FORMANTS {
            set_oral_formant_casc(&mut self.oral_formant_casc[i], f_parms, i)?;
        }

        // Adjust parallel branch:
        self.f_state.parallel_voicing_lin = db_to_lin(f_parms.parallel_voicing_db);
        self.f_state.parallel_aspiration_lin = db_to_lin(f_parms.parallel_aspiration_db);
        self.f_state.frication_lin = db_to_lin(f_parms.frication_db);
        self.f_state.parallel_bypass_lin = db_to_lin(f_parms.parallel_bypass_db);
        set_nasal_formant_par(&mut self.nasal_formant_par, f_parms)?;
        for i in 0..MAX_ORAL_FORMANTS {
            set_oral_formant_par(&mut self.oral_formant_par[i], m_parms, f_parms, i)?;
        }

        Ok(())
    }

    fn start_glottal_source_period(&mut self) -> KlResult<()> {
        self.glottal_source
            .start_period(self.p_state.as_ref().unwrap().open_phase_length)?;
        Ok(())
    }
}

fn set_tilt_filter(tilt_filter: &mut LpFilter1, tilt_db: f32) -> KlResult<()> {
    if tilt_db == 0. {
        tilt_filter.set_passthrough();
    } else {
        let tilt = db_to_lin(tilt_db);
        tilt_filter.set(3000., tilt, None)?;
    }
    Ok(())
}
fn set_nasal_formant_casc(
    nasal_formant_casc: &mut Resonator,
    f_parms: &FrameParms,
) -> KlResult<()> {
    if !valid_freq(f_parms.nasal_formant_freq) || !valid_freq(f_parms.nasal_formant_bw) {
        nasal_formant_casc.set_passthrough();
    } else {
        nasal_formant_casc.set(f_parms.nasal_formant_freq, f_parms.nasal_formant_bw, None)?;
    }
    Ok(())
}
fn set_nasal_antiformant_casc(
    nasal_antiformant_casc: &mut AntiResonator,
    f_parms: &FrameParms,
) -> KlResult<()> {
    if !valid_freq(f_parms.nasal_antiformant_freq) || !valid_freq(f_parms.nasal_antiformant_bw) {
        nasal_antiformant_casc.set_passthrough();
    } else {
        nasal_antiformant_casc.set(f_parms.nasal_antiformant_freq, f_parms.nasal_antiformant_bw)?;
    }
    Ok(())
}
fn set_oral_formant_casc(
    oral_formant_casc: &mut Resonator,
    f_parms: &FrameParms,
    i: usize,
) -> KlResult<()> {
    if i >= MAX_ORAL_FORMANTS {
        return KlError::err("oral formant index out of range");
    }
    let f = f_parms.oral_formant_freq[i];
    let bw = f_parms.oral_formant_bw[i];
    if !valid_freq(f) || !valid_freq(bw) {
        oral_formant_casc.set_passthrough();
    } else {
        oral_formant_casc.set(f, bw, None)?;
    }
    Ok(())
}
fn set_nasal_formant_par(nasal_formant_par: &mut Resonator, f_parms: &FrameParms) -> KlResult<()> {
    let lin = db_to_lin(f_parms.nasal_formant_db);
    if !valid_freq(f_parms.nasal_formant_freq) || !valid_freq(f_parms.nasal_formant_bw) || lin == 0.
    {
        nasal_formant_par.set_mute();
    } else {
        nasal_formant_par.set(f_parms.nasal_formant_freq, f_parms.nasal_formant_bw, None)?;
        nasal_formant_par.adjust_peak_gain(lin)?;
    }
    Ok(())
}
fn set_oral_formant_par(
    oral_formant_par: &mut Resonator,
    m_parms: &MainParms,
    f_parms: &FrameParms,
    i: usize,
) -> KlResult<()> {
    let formant = i + 1;
    let (f, bw, db) = if i < MAX_ORAL_FORMANTS {
        (
            f_parms.oral_formant_freq[i],
            f_parms.oral_formant_bw[i],
            f_parms.oral_formant_db[i],
        )
    } else {
        (f32::NAN, f32::NAN, f32::NAN)
    };
    let peak_gain = db_to_lin(db);
    // Klatt used the following linear factors to adjust the levels of the parallel formant
    // resonators so that they have a similar effect as the cascade versions:
    //   F1: 0.4, F2: 0.15, F3: 0.06, F4: 0.04, F5: 0.022, F6: 0.03, Nasal: 0.6
    // We are not doing this here, because then the output of the parallel branch would no longer
    // match the specified formant levels. Instead, we use the specified dB value to set the peak gain
    // instead of taking it as the DC gain.
    if !valid_freq(f) || !valid_freq(bw) || peak_gain == 0. {
        oral_formant_par.set_mute();
    } else {
        oral_formant_par.set(f, bw, None)?;
        let w = 2. * PI * f / m_parms.sample_rate as f32;
        // gain of differencing filter
        let diff_gain = (2. - 2. * w.cos()).sqrt();
        // compensate differencing filter for F2 to F6
        let filter_gain = if formant >= 2 {
            peak_gain / diff_gain
        } else {
            peak_gain
        };
        oral_formant_par.adjust_peak_gain(filter_gain)?;
    }
    Ok(())
}
fn adjust_signal_gain(buf: &mut [f32], target_rms: f32) {
    let rms = compute_rms(buf);
    if rms == 0. {
        return;
    }
    let r = target_rms / rms;
    for x in buf {
        *x *= r;
    }
}
fn compute_rms(buf: &[f32]) -> f32 {
    (buf.iter().map(|&x| x.powi(2)).sum::<f32>() / buf.len() as f32).sqrt()
}

/// Generates a sound that consists of multiple frames.
pub fn generate_sound<R: Rng + Clone + 'static>(
    m_parms: &MainParms,
    f_parms_a: &[FrameParms],
    rng: R,
) -> KlResult<Vec<f32>> {
    let mut generator = Generator::new(m_parms, rng)?;
    let out_buf_len = f_parms_a
        .iter()
        .map(|f_parms| (f_parms.duration * m_parms.sample_rate as f32).round() as usize)
        .sum();
    let mut out_buf = vec![0.; out_buf_len];
    let mut out_buf_pos = 0;
    for f_parms in f_parms_a {
        let frame_len = (f_parms.duration * m_parms.sample_rate as f32).round() as usize;
        let frame_buf = &mut out_buf[out_buf_pos..out_buf_pos + frame_len];
        generator.generate_frame(f_parms, frame_buf)?;
        out_buf_pos += frame_len;
    }
    Ok(out_buf)
}

// Transfer function

const EPS: f32 = 1e-10;

/// Returns the polynomial coefficients of the overall filter transfer function in the z-plane.
/// The returned array contains the top and bottom coefficients of the rational fraction, ordered in ascending powers.
pub fn get_vocal_tract_transfer_function_coefficients(
    m_parms: &MainParms,
    f_parms: &FrameParms,
) -> KlResult<Fraction> {
    let mut voice = fraction!([1], [1]); // glottal source

    let mut tilt_filter = LpFilter1::new(m_parms.sample_rate);
    set_tilt_filter(&mut tilt_filter, f_parms.tilt_db)?;
    let tilt_trans = tilt_filter.get_transfer_function_coefficients();
    voice = multiply_fractions(&voice, &tilt_trans, Some(EPS))?;

    let cascade_trans = if f_parms.cascade_enabled {
        get_cascade_branch_transfer_function_coefficients(m_parms, f_parms)?
    } else {
        fraction!([0], [1])
    };
    let parallel_trans = if f_parms.parallel_enabled {
        get_parallel_branch_transfer_function_coefficients(m_parms, f_parms)?
    } else {
        fraction!([0], [1])
    };
    let branches_trans = add_fractions(&cascade_trans, &parallel_trans, Some(EPS))?;
    voice = multiply_fractions(&voice, &branches_trans, Some(EPS))?;

    let mut output_lp_filter = Resonator::new(m_parms.sample_rate);
    output_lp_filter.set(0., m_parms.sample_rate as f32 / 2., None)?;
    let output_lp_trans = output_lp_filter.get_transfer_function_coefficients();
    voice = multiply_fractions(&voice, &output_lp_trans, Some(EPS))?;

    let gain_lin = db_to_lin(if f_parms.gain_db.is_nan() {
        0.
    } else {
        f_parms.gain_db
    });
    voice = multiply_fractions(&voice, &fraction!([gain_lin], [1]), Some(EPS))?;

    Ok(voice)
}

fn get_cascade_branch_transfer_function_coefficients(
    m_parms: &MainParms,
    f_parms: &FrameParms,
) -> KlResult<Fraction> {
    let cascade_voicing_lin = db_to_lin(f_parms.cascade_voicing_db);
    let mut v = fraction!([cascade_voicing_lin], [1]);

    let mut nasal_antiformant_casc = AntiResonator::new(m_parms.sample_rate);
    set_nasal_antiformant_casc(&mut nasal_antiformant_casc, f_parms)?;
    let nasal_antiformant_trans = nasal_antiformant_casc.get_transfer_function_coefficients();
    v = multiply_fractions(&v, &nasal_antiformant_trans, Some(EPS))?;

    let mut nasal_formant_casc = Resonator::new(m_parms.sample_rate);
    set_nasal_formant_casc(&mut nasal_formant_casc, f_parms)?;
    let nasal_formant_trans = nasal_formant_casc.get_transfer_function_coefficients();
    v = multiply_fractions(&v, &nasal_formant_trans, Some(EPS))?;

    for i in 0..MAX_ORAL_FORMANTS {
        let mut oral_formant_casc = Resonator::new(m_parms.sample_rate);
        set_oral_formant_casc(&mut oral_formant_casc, f_parms, i)?;
        let oral_formant_casc_trans = oral_formant_casc.get_transfer_function_coefficients();
        v = multiply_fractions(&v, &oral_formant_casc_trans, Some(EPS))?;
    }

    Ok(v)
}

fn get_parallel_branch_transfer_function_coefficients(
    m_parms: &MainParms,
    f_parms: &FrameParms,
) -> KlResult<Fraction> {
    let parallel_voicing_lin = db_to_lin(f_parms.parallel_voicing_db);
    let source = fraction!([parallel_voicing_lin], [1]);

    let differencing_filter_par = DifferencingFilter::new();
    let differencing_filter_trans = differencing_filter_par.get_transfer_function_coefficients();
    let source2 = multiply_fractions(&source, &differencing_filter_trans, Some(EPS))?;

    let mut v = fraction!([0], [1]);

    let mut nasal_formant_par = Resonator::new(m_parms.sample_rate);
    set_nasal_formant_par(&mut nasal_formant_par, f_parms)?;
    let nasal_formant_trans = nasal_formant_par.get_transfer_function_coefficients();
    v = add_fractions(
        &v,
        &multiply_fractions(&source, &nasal_formant_trans, Some(EPS))?,
        Some(EPS),
    )?;

    for i in 0..MAX_ORAL_FORMANTS {
        let mut oral_formant_par = Resonator::new(m_parms.sample_rate);
        set_oral_formant_par(&mut oral_formant_par, m_parms, f_parms, i)?;
        let oral_formant_trans = oral_formant_par.get_transfer_function_coefficients();
        let formant_in = if i == 0 { &source } else { &source2 }; // F1 is applied to source, F2 to F6 are applied to difference
        let formant_out = multiply_fractions(formant_in, &oral_formant_trans, Some(EPS))?;
        let alternating_sign = if i % 2 == 0 { 1. } else { -1. };
        let v2 = multiply_fractions(&formant_out, &fraction!([alternating_sign], [1]), Some(EPS))?;
        v = add_fractions(&v, &v2, Some(EPS))?;
    }

    let parallel_bypass_lin = db_to_lin(f_parms.parallel_bypass_db);
    let parallel_bypass =
        multiply_fractions(&source2, &fraction!([parallel_bypass_lin], [1]), Some(EPS))?;
    v = add_fractions(&v, &parallel_bypass, Some(EPS))?;

    Ok(v)
}
