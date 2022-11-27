use super::Matrix;
use std::collections::HashMap;

lazy_static! {
    pub static ref MATRIX: Matrix = HashMap::from([
        (
            65,
            HashMap::from([
                (65, 1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            82,
            HashMap::from([
                (65, -1),
                (82, 1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            78,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, 1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            66,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, 1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            68,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, 1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            67,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, 1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            81,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, 1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            90,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, 1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            69,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, 1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            71,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, 1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            72,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, 1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            73,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, 1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            76,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, 1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            75,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, 1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            77,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, 1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            70,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, 1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            80,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, 1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            83,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, 1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            84,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, 1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            87,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, 1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            89,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, 1),
                (86, -1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            86,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, 1),
                (88, -1),
                (42, -1),
            ])
        ),
        (
            88,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, 0),
                (42, -1),
            ])
        ),
        (
            42,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (66, -1),
                (68, -1),
                (67, -1),
                (81, -1),
                (90, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -1),
                (89, -1),
                (86, -1),
                (88, -1),
                (42, 0),
            ])
        ),
    ]);
}