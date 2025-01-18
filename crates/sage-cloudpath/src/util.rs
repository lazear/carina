use crate::{read_and_execute, tdf::BrukerSpectrumProcessor, Error};
use sage_core::spectrum::RawSpectrum;
use serde::Serialize;
use tokio::io::AsyncReadExt;

#[derive(Debug, PartialEq, Eq)]
enum FileFormat {
    MzML,
    MGF,
    TDF,
    Unidentified,
}

const BRUKER_EXTENSIONS: [&str; 5] = [".d", ".tdf", ".tdf_bin", "ms2", "raw"];

fn is_bruker(path: &str) -> bool {
    BRUKER_EXTENSIONS.iter().any(|ext| {
        if path.ends_with(std::path::MAIN_SEPARATOR) {
            path.strip_suffix(std::path::MAIN_SEPARATOR)
                .unwrap()
                .ends_with(ext)
        } else {
            path.ends_with(ext)
        }
    })
}

fn identify_format(s: &str) -> FileFormat {
    let path_lower = s.to_lowercase();
    if path_lower.ends_with(".mgf.gz") || path_lower.ends_with(".mgf") {
        FileFormat::MGF
    } else if is_bruker(&path_lower) {
        FileFormat::TDF
    } else if path_lower.ends_with(".mzml.gz") || path_lower.ends_with(".mzml") {
        FileFormat::MzML
    } else {
        FileFormat::Unidentified
    }
}

pub fn read_spectra<S: AsRef<str>>(
    path: S,
    file_id: usize,
    sn: Option<u8>,
    bruker_processor: BrukerSpectrumProcessor,
    requires_ms1: bool,
) -> Result<Vec<RawSpectrum>, Error> {
    match identify_format(path.as_ref()) {
        FileFormat::MzML => read_mzml(path, file_id, sn),
        FileFormat::MGF => read_mgf(path, file_id),
        FileFormat::TDF => read_tdf(path, file_id, bruker_processor, requires_ms1),
        FileFormat::Unidentified => panic!("Unable to get type for '{}'", path.as_ref()), // read_mzml(path, file_id, sn),
    }
}

pub fn read_mzml<S: AsRef<str>>(
    s: S,
    file_id: usize,
    signal_to_noise: Option<u8>,
) -> Result<Vec<RawSpectrum>, Error> {
    read_and_execute(s, |bf| async move {
        Ok(crate::mzml::MzMLReader::with_file_id(file_id)
            .set_signal_to_noise(signal_to_noise)
            .parse(bf)
            .await?)
    })
}

pub fn read_tdf<S: AsRef<str>>(
    s: S,
    file_id: usize,
    bruker_spectrum_processor: BrukerSpectrumProcessor,
    requires_ms1: bool,
) -> Result<Vec<RawSpectrum>, Error> {
    let res = crate::tdf::TdfReader.parse(s, file_id, bruker_spectrum_processor, requires_ms1);
    match res {
        Ok(t) => Ok(t),
        Err(e) => Err(Error::TDF(e)),
    }
}

pub fn read_mgf<S: AsRef<str>>(path: S, file_id: usize) -> Result<Vec<RawSpectrum>, Error> {
    read_and_execute(path, |mut bf| async move {
        let mut contents = String::new();
        bf.read_to_string(&mut contents)
            .await
            .map_err(crate::Error::IO)?;
        let res = crate::mgf::MgfReader::with_file_id(file_id).parse(contents);
        match res {
            Ok(m) => Ok(m),
            Err(e) => Err(Error::MGF(e)),
        }
    })
}

pub fn read_fasta<S>(
    path: S,
    decoy_tag: S,
    generate_decoys: bool,
) -> Result<sage_core::fasta::Fasta, Error>
where
    S: AsRef<str>,
{
    read_and_execute(path, |mut bf| async move {
        let mut contents = String::new();
        bf.read_to_string(&mut contents)
            .await
            .map_err(crate::Error::IO)?;
        Ok(sage_core::fasta::Fasta::parse(
            contents,
            decoy_tag.as_ref(),
            generate_decoys,
        ))
    })
}

pub fn read_json<S, T>(path: S) -> Result<T, Error>
where
    S: AsRef<str>,
    T: for<'de> serde::Deserialize<'de>,
{
    read_and_execute(path, |mut bf| async move {
        let mut contents = String::new();
        bf.read_to_string(&mut contents).await?;
        Ok(serde_json::from_str(&contents)?)
    })
}

/// Send telemetry data
pub fn send_data<T>(url: &str, data: &T) -> Result<(), Box<dyn std::error::Error + Send + Sync>>
where
    T: Serialize,
{
    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()?;

    rt.block_on(async {
        let client = reqwest::ClientBuilder::default().https_only(true).build()?;
        let res = client.post(url).json(data).send().await?;
        res.error_for_status()?;
        Ok(())
    })
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_identify_format() {
        assert_eq!(identify_format("foo.mzml"), FileFormat::MzML);
        assert_eq!(identify_format("foo.mzML"), FileFormat::MzML);
        assert_eq!(identify_format("foo.mgf"), FileFormat::MGF);
        assert_eq!(identify_format("foo.mgf.gz"), FileFormat::MGF);
        assert_eq!(identify_format("foo.tdf"), FileFormat::TDF);
        assert_eq!(identify_format("./tomato/foo.d"), FileFormat::TDF);
        assert_eq!(identify_format("./tomato/foo.d/"), FileFormat::TDF);
    }
}
