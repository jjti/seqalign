use super::Matrix;
use std::collections::HashMap;

lazy_static! {
    pub static ref MATRIX: Matrix = HashMap::from([
        (
            65,
            HashMap::from([
                (65, 3),
                (82, -3),
                (78, -1),
                (68, 0),
                (67, -3),
                (81, -1),
                (69, 0),
                (71, 1),
                (72, -3),
                (73, -1),
                (76, -3),
                (75, -2),
                (77, -2),
                (70, -4),
                (80, 1),
                (83, 1),
                (84, 1),
                (87, -7),
                (89, -4),
                (86, 0),
                (66, 0),
                (90, -1),
                (88, -1),
                (42, -8),
            ])
        ),
        (
            82,
            HashMap::from([
                (65, -3),
                (82, 6),
                (78, -1),
                (68, -3),
                (67, -4),
                (81, 1),
                (69, -3),
                (71, -4),
                (72, 1),
                (73, -2),
                (76, -4),
                (75, 2),
                (77, -1),
                (70, -5),
                (80, -1),
                (83, -1),
                (84, -2),
                (87, 1),
                (89, -5),
                (86, -3),
                (66, -2),
                (90, -1),
                (88, -2),
                (42, -8),
            ])
        ),
        (
            78,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, 4),
                (68, 2),
                (67, -5),
                (81, 0),
                (69, 1),
                (71, 0),
                (72, 2),
                (73, -2),
                (76, -4),
                (75, 1),
                (77, -3),
                (70, -4),
                (80, -2),
                (83, 1),
                (84, 0),
                (87, -4),
                (89, -2),
                (86, -3),
                (66, 3),
                (90, 0),
                (88, -1),
                (42, -8),
            ])
        ),
        (
            68,
            HashMap::from([
                (65, 0),
                (82, -3),
                (78, 2),
                (68, 5),
                (67, -7),
                (81, 1),
                (69, 3),
                (71, 0),
                (72, 0),
                (73, -3),
                (76, -5),
                (75, -1),
                (77, -4),
                (70, -7),
                (80, -3),
                (83, 0),
                (84, -1),
                (87, -8),
                (89, -5),
                (86, -3),
                (66, 4),
                (90, 3),
                (88, -2),
                (42, -8),
            ])
        ),
        (
            67,
            HashMap::from([
                (65, -3),
                (82, -4),
                (78, -5),
                (68, -7),
                (67, 9),
                (81, -7),
                (69, -7),
                (71, -4),
                (72, -4),
                (73, -3),
                (76, -7),
                (75, -7),
                (77, -6),
                (70, -6),
                (80, -4),
                (83, 0),
                (84, -3),
                (87, -8),
                (89, -1),
                (86, -3),
                (66, -6),
                (90, -7),
                (88, -4),
                (42, -8),
            ])
        ),
        (
            81,
            HashMap::from([
                (65, -1),
                (82, 1),
                (78, 0),
                (68, 1),
                (67, -7),
                (81, 6),
                (69, 2),
                (71, -3),
                (72, 3),
                (73, -3),
                (76, -2),
                (75, 0),
                (77, -1),
                (70, -6),
                (80, 0),
                (83, -2),
                (84, -2),
                (87, -6),
                (89, -5),
                (86, -3),
                (66, 0),
                (90, 4),
                (88, -1),
                (42, -8),
            ])
        ),
        (
            69,
            HashMap::from([
                (65, 0),
                (82, -3),
                (78, 1),
                (68, 3),
                (67, -7),
                (81, 2),
                (69, 5),
                (71, -1),
                (72, -1),
                (73, -3),
                (76, -4),
                (75, -1),
                (77, -3),
                (70, -7),
                (80, -2),
                (83, -1),
                (84, -2),
                (87, -8),
                (89, -5),
                (86, -3),
                (66, 3),
                (90, 4),
                (88, -1),
                (42, -8),
            ])
        ),
        (
            71,
            HashMap::from([
                (65, 1),
                (82, -4),
                (78, 0),
                (68, 0),
                (67, -4),
                (81, -3),
                (69, -1),
                (71, 5),
                (72, -4),
                (73, -4),
                (76, -5),
                (75, -3),
                (77, -4),
                (70, -5),
                (80, -2),
                (83, 1),
                (84, -1),
                (87, -8),
                (89, -6),
                (86, -2),
                (66, 0),
                (90, -2),
                (88, -2),
                (42, -8),
            ])
        ),
        (
            72,
            HashMap::from([
                (65, -3),
                (82, 1),
                (78, 2),
                (68, 0),
                (67, -4),
                (81, 3),
                (69, -1),
                (71, -4),
                (72, 7),
                (73, -4),
                (76, -3),
                (75, -2),
                (77, -4),
                (70, -3),
                (80, -1),
                (83, -2),
                (84, -3),
                (87, -3),
                (89, -1),
                (86, -3),
                (66, 1),
                (90, 1),
                (88, -2),
                (42, -8),
            ])
        ),
        (
            73,
            HashMap::from([
                (65, -1),
                (82, -2),
                (78, -2),
                (68, -3),
                (67, -3),
                (81, -3),
                (69, -3),
                (71, -4),
                (72, -4),
                (73, 6),
                (76, 1),
                (75, -3),
                (77, 1),
                (70, 0),
                (80, -3),
                (83, -2),
                (84, 0),
                (87, -6),
                (89, -2),
                (86, 3),
                (66, -3),
                (90, -3),
                (88, -1),
                (42, -8),
            ])
        ),
        (
            76,
            HashMap::from([
                (65, -3),
                (82, -4),
                (78, -4),
                (68, -5),
                (67, -7),
                (81, -2),
                (69, -4),
                (71, -5),
                (72, -3),
                (73, 1),
                (76, 5),
                (75, -4),
                (77, 3),
                (70, 0),
                (80, -3),
                (83, -4),
                (84, -3),
                (87, -3),
                (89, -2),
                (86, 1),
                (66, -4),
                (90, -3),
                (88, -2),
                (42, -8),
            ])
        ),
        (
            75,
            HashMap::from([
                (65, -2),
                (82, 2),
                (78, 1),
                (68, -1),
                (67, -7),
                (81, 0),
                (69, -1),
                (71, -3),
                (72, -2),
                (73, -3),
                (76, -4),
                (75, 5),
                (77, 0),
                (70, -7),
                (80, -2),
                (83, -1),
                (84, -1),
                (87, -5),
                (89, -5),
                (86, -4),
                (66, 0),
                (90, -1),
                (88, -2),
                (42, -8),
            ])
        ),
        (
            77,
            HashMap::from([
                (65, -2),
                (82, -1),
                (78, -3),
                (68, -4),
                (67, -6),
                (81, -1),
                (69, -3),
                (71, -4),
                (72, -4),
                (73, 1),
                (76, 3),
                (75, 0),
                (77, 8),
                (70, -1),
                (80, -3),
                (83, -2),
                (84, -1),
                (87, -6),
                (89, -4),
                (86, 1),
                (66, -4),
                (90, -2),
                (88, -2),
                (42, -8),
            ])
        ),
        (
            70,
            HashMap::from([
                (65, -4),
                (82, -5),
                (78, -4),
                (68, -7),
                (67, -6),
                (81, -6),
                (69, -7),
                (71, -5),
                (72, -3),
                (73, 0),
                (76, 0),
                (75, -7),
                (77, -1),
                (70, 8),
                (80, -5),
                (83, -3),
                (84, -4),
                (87, -1),
                (89, 4),
                (86, -3),
                (66, -5),
                (90, -6),
                (88, -3),
                (42, -8),
            ])
        ),
        (
            80,
            HashMap::from([
                (65, 1),
                (82, -1),
                (78, -2),
                (68, -3),
                (67, -4),
                (81, 0),
                (69, -2),
                (71, -2),
                (72, -1),
                (73, -3),
                (76, -3),
                (75, -2),
                (77, -3),
                (70, -5),
                (80, 6),
                (83, 1),
                (84, -1),
                (87, -7),
                (89, -6),
                (86, -2),
                (66, -2),
                (90, -1),
                (88, -2),
                (42, -8),
            ])
        ),
        (
            83,
            HashMap::from([
                (65, 1),
                (82, -1),
                (78, 1),
                (68, 0),
                (67, 0),
                (81, -2),
                (69, -1),
                (71, 1),
                (72, -2),
                (73, -2),
                (76, -4),
                (75, -1),
                (77, -2),
                (70, -3),
                (80, 1),
                (83, 3),
                (84, 2),
                (87, -2),
                (89, -3),
                (86, -2),
                (66, 0),
                (90, -1),
                (88, -1),
                (42, -8),
            ])
        ),
        (
            84,
            HashMap::from([
                (65, 1),
                (82, -2),
                (78, 0),
                (68, -1),
                (67, -3),
                (81, -2),
                (69, -2),
                (71, -1),
                (72, -3),
                (73, 0),
                (76, -3),
                (75, -1),
                (77, -1),
                (70, -4),
                (80, -1),
                (83, 2),
                (84, 4),
                (87, -6),
                (89, -3),
                (86, 0),
                (66, 0),
                (90, -2),
                (88, -1),
                (42, -8),
            ])
        ),
        (
            87,
            HashMap::from([
                (65, -7),
                (82, 1),
                (78, -4),
                (68, -8),
                (67, -8),
                (81, -6),
                (69, -8),
                (71, -8),
                (72, -3),
                (73, -6),
                (76, -3),
                (75, -5),
                (77, -6),
                (70, -1),
                (80, -7),
                (83, -2),
                (84, -6),
                (87, 12),
                (89, -2),
                (86, -8),
                (66, -6),
                (90, -7),
                (88, -5),
                (42, -8),
            ])
        ),
        (
            89,
            HashMap::from([
                (65, -4),
                (82, -5),
                (78, -2),
                (68, -5),
                (67, -1),
                (81, -5),
                (69, -5),
                (71, -6),
                (72, -1),
                (73, -2),
                (76, -2),
                (75, -5),
                (77, -4),
                (70, 4),
                (80, -6),
                (83, -3),
                (84, -3),
                (87, -2),
                (89, 8),
                (86, -3),
                (66, -3),
                (90, -5),
                (88, -3),
                (42, -8),
            ])
        ),
        (
            86,
            HashMap::from([
                (65, 0),
                (82, -3),
                (78, -3),
                (68, -3),
                (67, -3),
                (81, -3),
                (69, -3),
                (71, -2),
                (72, -3),
                (73, 3),
                (76, 1),
                (75, -4),
                (77, 1),
                (70, -3),
                (80, -2),
                (83, -2),
                (84, 0),
                (87, -8),
                (89, -3),
                (86, 5),
                (66, -3),
                (90, -3),
                (88, -1),
                (42, -8),
            ])
        ),
        (
            66,
            HashMap::from([
                (65, 0),
                (82, -2),
                (78, 3),
                (68, 4),
                (67, -6),
                (81, 0),
                (69, 3),
                (71, 0),
                (72, 1),
                (73, -3),
                (76, -4),
                (75, 0),
                (77, -4),
                (70, -5),
                (80, -2),
                (83, 0),
                (84, 0),
                (87, -6),
                (89, -3),
                (86, -3),
                (66, 4),
                (90, 2),
                (88, -1),
                (42, -8),
            ])
        ),
        (
            90,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, 0),
                (68, 3),
                (67, -7),
                (81, 4),
                (69, 4),
                (71, -2),
                (72, 1),
                (73, -3),
                (76, -3),
                (75, -1),
                (77, -2),
                (70, -6),
                (80, -1),
                (83, -1),
                (84, -2),
                (87, -7),
                (89, -5),
                (86, -3),
                (66, 2),
                (90, 4),
                (88, -1),
                (42, -8),
            ])
        ),
        (
            88,
            HashMap::from([
                (65, -1),
                (82, -2),
                (78, -1),
                (68, -2),
                (67, -4),
                (81, -1),
                (69, -1),
                (71, -2),
                (72, -2),
                (73, -1),
                (76, -2),
                (75, -2),
                (77, -2),
                (70, -3),
                (80, -2),
                (83, -1),
                (84, -1),
                (87, -5),
                (89, -3),
                (86, -1),
                (66, -1),
                (90, -1),
                (88, -2),
                (42, -8),
            ])
        ),
        (
            42,
            HashMap::from([
                (65, -8),
                (82, -8),
                (78, -8),
                (68, -8),
                (67, -8),
                (81, -8),
                (69, -8),
                (71, -8),
                (72, -8),
                (73, -8),
                (76, -8),
                (75, -8),
                (77, -8),
                (70, -8),
                (80, -8),
                (83, -8),
                (84, -8),
                (87, -8),
                (89, -8),
                (86, -8),
                (66, -8),
                (90, -8),
                (88, -8),
                (42, 1),
            ])
        ),
    ]);
}
