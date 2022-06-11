use clap::Parser;
use reqwest::StatusCode;
use rust_htslib::bcf::{Read, Reader, Record};
use rust_htslib::errors::Result;
use serde::{Deserialize, Serialize};
use std::str;

const GENE_SYMBOL: usize = 3;
const TYPE: usize = 5;
const HGVS: usize = 9;
const NS: &'static str = "NS";
const DP: &'static str = "DP";
const GQ: &'static str = "GQ";

#[derive(Serialize, Deserialize, Debug)]
struct Distribution {
    min: f32,
    q25: f32,
    median: f32,
    q75: f32,
    max: f32,
    mean: f32,
}

#[derive(Serialize, Deserialize, Debug)]
struct Variant {
    id: Option<String>,
    #[serde(rename = "datasetId")]
    dataset_id: String,
    #[serde(rename = "totalSamples")]
    total_samples: u32,
    #[serde(rename = "assemblyId")]
    assembly_id: String,
    #[serde(rename = "snpIds")]
    snp_ids: Option<Vec<String>>,
    #[serde(rename = "referenceName")]
    reference_name: String,
    start: i64,
    #[serde(rename = "referenceBases")]
    reference_bases: String,
    #[serde(rename = "alternateBases")]
    alternate_bases: Vec<String>,
    #[serde(rename = "geneSymbol")]
    gene_symbol: Option<Vec<String>>,
    #[serde(rename = "alleleFrequency")]
    allele_frequency: Vec<f32>,
    #[serde(rename = "sampleCount")]
    sample_count: Option<i32>,
    coverage: Distribution,
    #[serde(rename = "genotypeQuality")]
    genotype_quality: Distribution,
    clnsig: Option<String>,
    hgvs: Option<Vec<String>>,
    #[serde(rename = "type")]
    variant_type: Option<Vec<String>>,
}

fn get_snp_ids(record: &Record) -> Option<Vec<String>> {
    let id = record.id();
    let id = str::from_utf8(&id).unwrap();
    if id == "." {
        return None;
    }
    Some(id.split(';').map(|x| x.to_string()).collect())
}

fn get_allele_frequency(record: &Record) -> Result<Option<Vec<f32>>> {
    Ok(record.info("AF".as_bytes()).float()?.map(|x| x.to_vec()))
}

#[derive(Parser)]
#[clap(version = "0.1.0", author = "Welliton de Souza <well309@gmail.com>")]
struct Opts {
    #[clap(
        long,
        default_value = "http://localhost:8080",
        help = "URL to BraVE server"
    )]
    host: String,
    #[clap(long, help = "Dataset name")]
    dataset: String,
    #[clap(long, help = "Genome assembly version")]
    assembly: String,
    #[clap(long, default_value = "admin", help = "User name")]
    username: String,
    #[clap(long, help = "Password")]
    password: Option<String>,
    #[clap(long, help = "Don't filter variants by FILTER column")]
    dont_filter: bool,
    #[clap(long, help = "Just check VCF without connecting to server")]
    dryrun: bool,
    #[clap(long, help = "Print variant data to stderr")]
    debug: bool,
    #[clap(long, help = "Disable SSL certification verification")]
    disable_ssl: bool,
    vcf_file: String,
}

fn main() {
    let opts: Opts = Opts::parse();

    let dataset_id = opts.dataset;
    let assemble_id = opts.assembly;
    let do_filter = !opts.dont_filter;
    let path = opts.vcf_file;
    let host = opts.host;
    let username = opts.username;
    let password = opts.password;
    let dryrun = opts.dryrun;
    let debug = opts.debug;
    let disable_ssl = opts.disable_ssl;

    let mut bcf = Reader::from_path(path).expect("Error opening file.");

    let total_samples = bcf.header().sample_count();

    let has_ns = bcf.header().info_type(NS.as_bytes()).is_ok();

    let client = reqwest::blocking::Client::builder().danger_accept_invalid_certs(!disable_ssl);

    let url = format!("{}/variants", host);

    let mut total_variants: u32 = 0;
    let mut passed_variants: u32 = 0;

    for record in bcf.records() {
        let record = record.unwrap();

        total_variants += 1;

        if do_filter && !record.has_filter("PASS".as_bytes()) {
            continue;
        }

        passed_variants += 1;

        let snp_ids = get_snp_ids(&record);
        let allele_frequency: Vec<f32> = get_allele_frequency(&record).unwrap().unwrap_or_default();
        let coverage = calc_distribution(&record, DP);
        let genotype_quality = calc_distribution(&record, GQ);
        let start = record.pos() + 1;

        let rid = record
            .rid()
            .unwrap_or_else(|| panic!("Missing CHROM at position {}", start));

        let reference_name = record
            .header()
            .rid2name(rid)
            .map(|x| str::from_utf8(x).unwrap().to_string())
            .unwrap();

        let reference_bases = record
            .alleles()
            .first()
            .map(|x| str::from_utf8(x).unwrap().to_string())
            .unwrap_or_else(|| panic!("Missing REF at position {}", start));

        let alternate_bases = record
            .alleles()
            .iter()
            .skip(1)
            .map(|x| str::from_utf8(x).unwrap().to_string())
            .collect();

        let clnsig = get_info_field(&record, "CLNSIG").map(|x| x[0].to_owned());

        let sample_count = if has_ns {
            record.info(NS.as_bytes()).integer().unwrap().map(|x| x[0])
        } else {
            None
        };

        let maybe_ann = get_info_field(&record, "ANN");
        let (gene_symbol, variant_type, hgvs) = if let Some(ann) = maybe_ann {
            let fields: Vec<Vec<String>> = ann.iter().map(|x| split_ann(x)).collect();
            let gene_symbol = get_field(&fields, GENE_SYMBOL);
            let variant_type = get_field(&fields, TYPE);
            let hgvs = get_field(&fields, HGVS);
            (Some(gene_symbol), Some(variant_type), Some(hgvs))
        } else {
            (None, None, None)
        };

        let v = Variant {
            id: None,
            dataset_id: dataset_id.to_string(),
            total_samples,
            assembly_id: assemble_id.to_string(),
            snp_ids,
            reference_name,
            start,
            reference_bases,
            alternate_bases,
            gene_symbol,
            allele_frequency,
            sample_count,
            coverage,
            genotype_quality,
            clnsig,
            hgvs,
            variant_type,
        };

        if debug {
            eprintln!("{:?}", v);
        }

        if dryrun {
            continue;
        }

        let res = client
            .post(&url)
            .basic_auth(&username, password.as_ref())
            .json(&v)
            .send()
            .unwrap();
        assert_eq!(res.status(), StatusCode::CREATED, "{}", res.text().unwrap());
    }

    println!("Total variants: {}", total_variants);
    if do_filter {
        println!("Passed variants: {}", passed_variants);
    }
}

fn calc_distribution(record: &Record, tag: &str) -> Distribution {
    let mut values: Vec<i64> = record
        .format(tag.as_bytes())
        .integer()
        .unwrap()
        .iter()
        .map(|x| x[0] as i64)
        .collect();

    values.sort();

    Distribution {
        min: values[0] as f32,
        q25: percentile(&values, 25_f32),
        median: percentile(&values, 50_f32),
        q75: percentile(&values, 75_f32),
        max: values[values.len() - 1] as f32,
        mean: values.iter().sum::<i64>() as f32 / values.len() as f32,
    }
}

fn split_ann(ann: &str) -> Vec<String> {
    ann.split('|').map(|field| field.to_string()).collect()
}

fn get_field(fields: &Vec<Vec<String>>, index: usize) -> Vec<String> {
    fields
        .iter()
        .map(|x| x[index].to_string())
        .map(|x| x.to_string())
        .collect()
}

fn get_info_field(record: &Record, tag: &str) -> Option<Vec<String>> {
    let info = record.info(tag.as_bytes()).string().unwrap()?;
    Some(
        info.iter()
            .map(|y| str::from_utf8(y).unwrap().to_string())
            .collect(),
    )
}

// From https://doc.rust-lang.org/src/test/stats.rs.html
fn percentile(sorted_values: &[i64], pct: f32) -> f32 {
    assert!(!sorted_values.is_empty());

    if sorted_values.len() == 1 {
        return sorted_values[0] as f32;
    }

    let zero: f32 = 0.0;
    assert!(zero <= pct);

    let hundred = 100_f32;
    assert!(pct <= hundred);

    if pct == hundred {
        return sorted_values[sorted_values.len() - 1] as f32;
    }

    let length = (sorted_values.len() - 1) as f32;
    let rank = (pct / hundred) * length;
    let lrank = rank.floor();
    let d = rank - lrank;
    let n = lrank as usize;
    let lo = sorted_values[n] as f32;
    let hi = sorted_values[n + 1] as f32;
    lo + (hi - lo) * d
}
