use std::collections::HashMap;

/// Matrix is a single alignment matrix used in scoring an alignment.
///
/// Maps a char to another char and the corresponding substitution penalty.
pub type Matrix = HashMap<u8, HashMap<u8, i32>>;

/// All the modules beneath here were auto-generated.
#[allow(non_snake_case)]
pub mod BLOSUM100;
#[allow(non_snake_case)]
pub mod BLOSUM100_50;
#[allow(non_snake_case)]
pub mod BLOSUM30;
#[allow(non_snake_case)]
pub mod BLOSUM30_50;
#[allow(non_snake_case)]
pub mod BLOSUM35;
#[allow(non_snake_case)]
pub mod BLOSUM35_50;
#[allow(non_snake_case)]
pub mod BLOSUM40;
#[allow(non_snake_case)]
pub mod BLOSUM40_50;
#[allow(non_snake_case)]
pub mod BLOSUM45;
#[allow(non_snake_case)]
pub mod BLOSUM45_50;
#[allow(non_snake_case)]
pub mod BLOSUM50;
#[allow(non_snake_case)]
pub mod BLOSUM50_50;
#[allow(non_snake_case)]
pub mod BLOSUM55;
#[allow(non_snake_case)]
pub mod BLOSUM55_50;
#[allow(non_snake_case)]
pub mod BLOSUM60;
#[allow(non_snake_case)]
pub mod BLOSUM60_50;
#[allow(non_snake_case)]
pub mod BLOSUM62;
#[allow(non_snake_case)]
pub mod BLOSUM62_50;
#[allow(non_snake_case)]
pub mod BLOSUM65;
#[allow(non_snake_case)]
pub mod BLOSUM65_50;
#[allow(non_snake_case)]
pub mod BLOSUM70;
#[allow(non_snake_case)]
pub mod BLOSUM70_50;
#[allow(non_snake_case)]
pub mod BLOSUM75;
#[allow(non_snake_case)]
pub mod BLOSUM75_50;
#[allow(non_snake_case)]
pub mod BLOSUM80;
#[allow(non_snake_case)]
pub mod BLOSUM80_50;
#[allow(non_snake_case)]
pub mod BLOSUM85;
#[allow(non_snake_case)]
pub mod BLOSUM85_50;
#[allow(non_snake_case)]
pub mod BLOSUM90;
#[allow(non_snake_case)]
pub mod BLOSUM90_50;
#[allow(non_snake_case)]
pub mod BLOSUMN;
#[allow(non_snake_case)]
pub mod BLOSUMN_50;
#[allow(non_snake_case)]
pub mod DAYHOFF;
#[allow(non_snake_case)]
pub mod GONNET;
#[allow(non_snake_case)]
pub mod IDENTITY;
#[allow(non_snake_case)]
pub mod MATCH;
#[allow(non_snake_case)]
pub mod NUC_4_2;
#[allow(non_snake_case)]
pub mod NUC_4_4;
#[allow(non_snake_case)]
pub mod PAM10;
#[allow(non_snake_case)]
pub mod PAM100;
#[allow(non_snake_case)]
pub mod PAM110;
#[allow(non_snake_case)]
pub mod PAM120;
#[allow(non_snake_case)]
pub mod PAM120_CDI;
#[allow(non_snake_case)]
pub mod PAM130;
#[allow(non_snake_case)]
pub mod PAM140;
#[allow(non_snake_case)]
pub mod PAM150;
#[allow(non_snake_case)]
pub mod PAM160;
#[allow(non_snake_case)]
pub mod PAM160_CDI;
#[allow(non_snake_case)]
pub mod PAM170;
#[allow(non_snake_case)]
pub mod PAM180;
#[allow(non_snake_case)]
pub mod PAM190;
#[allow(non_snake_case)]
pub mod PAM20;
#[allow(non_snake_case)]
pub mod PAM200;
#[allow(non_snake_case)]
pub mod PAM200_CDI;
#[allow(non_snake_case)]
pub mod PAM210;
#[allow(non_snake_case)]
pub mod PAM220;
#[allow(non_snake_case)]
pub mod PAM230;
#[allow(non_snake_case)]
pub mod PAM240;
#[allow(non_snake_case)]
pub mod PAM250;
#[allow(non_snake_case)]
pub mod PAM250_CDI;
#[allow(non_snake_case)]
pub mod PAM260;
#[allow(non_snake_case)]
pub mod PAM270;
#[allow(non_snake_case)]
pub mod PAM280;
#[allow(non_snake_case)]
pub mod PAM290;
#[allow(non_snake_case)]
pub mod PAM30;
#[allow(non_snake_case)]
pub mod PAM300;
#[allow(non_snake_case)]
pub mod PAM310;
#[allow(non_snake_case)]
pub mod PAM320;
#[allow(non_snake_case)]
pub mod PAM330;
#[allow(non_snake_case)]
pub mod PAM340;
#[allow(non_snake_case)]
pub mod PAM350;
#[allow(non_snake_case)]
pub mod PAM360;
#[allow(non_snake_case)]
pub mod PAM370;
#[allow(non_snake_case)]
pub mod PAM380;
#[allow(non_snake_case)]
pub mod PAM390;
#[allow(non_snake_case)]
pub mod PAM40;
#[allow(non_snake_case)]
pub mod PAM400;
#[allow(non_snake_case)]
pub mod PAM40_CDI;
#[allow(non_snake_case)]
pub mod PAM410;
#[allow(non_snake_case)]
pub mod PAM420;
#[allow(non_snake_case)]
pub mod PAM430;
#[allow(non_snake_case)]
pub mod PAM440;
#[allow(non_snake_case)]
pub mod PAM450;
#[allow(non_snake_case)]
pub mod PAM460;
#[allow(non_snake_case)]
pub mod PAM470;
#[allow(non_snake_case)]
pub mod PAM480;
#[allow(non_snake_case)]
pub mod PAM490;
#[allow(non_snake_case)]
pub mod PAM50;
#[allow(non_snake_case)]
pub mod PAM500;
#[allow(non_snake_case)]
pub mod PAM60;
#[allow(non_snake_case)]
pub mod PAM70;
#[allow(non_snake_case)]
pub mod PAM80;
#[allow(non_snake_case)]
pub mod PAM80_CDI;
#[allow(non_snake_case)]
pub mod PAM90;
