use hound::{SampleFormat, WavSpec, WavWriter};
use klsyn::*;
use rand::rng;
use std::error::Error;

const DEMO_MAIN_PARMS: MainParms = MainParms {
    sample_rate: 44100,
    glottal_source_type: GlottalSourceType::Natural,
};
const DEMO_FRAME_PARMS: FrameParms = FrameParms {
    duration: 1.0,
    f0: 180.0,
    flutter_level: 0.25,
    open_phase_ratio: 0.7,
    breathiness_db: -25.0,
    tilt_db: 0.0,
    gain_db: f32::NAN,
    agc_rms_level: 0.2,
    nasal_formant_freq: f32::NAN,
    nasal_formant_bw: f32::NAN,
    oral_formant_freq: [520.0, 1006.0, 2831.0, 3168.0, 4135.0, 5020.0],
    oral_formant_bw: [76.0, 102.0, 72.0, 102.0, 816.0, 596.0],
    cascade_enabled: true,
    cascade_voicing_db: 0.0,
    cascade_aspiration_db: -25.0,
    cascade_aspiration_mod: 0.5,
    nasal_antiformant_freq: f32::NAN,
    nasal_antiformant_bw: f32::NAN,
    parallel_enabled: true,
    parallel_voicing_db: 0.0,
    parallel_aspiration_db: -25.0,
    parallel_aspiration_mod: 0.5,
    frication_db: -30.0,
    frication_mod: 0.5,
    parallel_bypass_db: -99.0,
    nasal_formant_db: f32::NAN,
    oral_formant_db: [0.0, -8.0, -15.0, -19.0, -30.0, -35.0],
};

fn main() -> Result<(), Box<dyn Error>> {
    let sound = generate_sound(&DEMO_MAIN_PARMS, &[DEMO_FRAME_PARMS], rng())?;

    let spec = WavSpec {
        channels: 1,
        sample_rate: DEMO_MAIN_PARMS.sample_rate as u32,
        bits_per_sample: 32,
        sample_format: SampleFormat::Float,
    };

    let mut writer = WavWriter::create("vowel.wav", spec)?;
    for sample in sound {
        writer.write_sample(sample)?;
    }

    println!(
        "Vocal tract transfer function coefficients: {}",
        get_vocal_tract_transfer_function_coefficients(&DEMO_MAIN_PARMS, &DEMO_FRAME_PARMS)?
    );

    Ok(())
}
