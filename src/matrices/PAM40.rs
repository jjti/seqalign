use super::Matrix;
use std::collections::HashMap;

lazy_static! {
    pub static ref MATRIX: Matrix = HashMap::from([
        (
            65,
            HashMap::from([
                (65, 6),
                (82, -6),
                (78, -3),
                (68, -3),
                (67, -6),
                (81, -3),
                (69, -2),
                (71, -1),
                (72, -6),
                (73, -4),
                (76, -5),
                (75, -6),
                (77, -4),
                (70, -7),
                (80, -1),
                (83, 0),
                (84, 0),
                (87, -12),
                (89, -7),
                (86, -2),
                (66, -3),
                (90, -2),
                (88, -3),
                (42, -15),
            ])
        ),
        (
            82,
            HashMap::from([
                (65, -6),
                (82, 8),
                (78, -5),
                (68, -9),
                (67, -7),
                (81, -1),
                (69, -8),
                (71, -8),
                (72, -1),
                (73, -5),
                (76, -8),
                (75, 1),
                (77, -3),
                (70, -8),
                (80, -3),
                (83, -2),
                (84, -5),
                (87, -1),
                (89, -9),
                (86, -7),
                (66, -6),
                (90, -3),
                (88, -5),
                (42, -15),
            ])
        ),
        (
            78,
            HashMap::from([
                (65, -3),
                (82, -5),
                (78, 7),
                (68, 2),
                (67, -9),
                (81, -3),
                (69, -1),
                (71, -2),
                (72, 1),
                (73, -4),
                (76, -6),
                (75, 0),
                (77, -7),
                (70, -8),
                (80, -5),
                (83, 0),
                (84, -1),
                (87, -7),
                (89, -4),
                (86, -7),
                (66, 6),
                (90, -2),
                (88, -3),
                (42, -15),
            ])
        ),
        (
            68,
            HashMap::from([
                (65, -3),
                (82, -9),
                (78, 2),
                (68, 7),
                (67, -12),
                (81, -2),
                (69, 3),
                (71, -3),
                (72, -3),
                (73, -6),
                (76, -11),
                (75, -4),
                (77, -9),
                (70, -13),
                (80, -7),
                (83, -3),
                (84, -4),
                (87, -13),
                (89, -10),
                (86, -7),
                (66, 6),
                (90, 2),
                (88, -5),
                (42, -15),
            ])
        ),
        (
            67,
            HashMap::from([
                (65, -6),
                (82, -7),
                (78, -9),
                (68, -12),
                (67, 9),
                (81, -12),
                (69, -12),
                (71, -8),
                (72, -7),
                (73, -5),
                (76, -13),
                (75, -12),
                (77, -12),
                (70, -11),
                (80, -7),
                (83, -2),
                (84, -7),
                (87, -14),
                (89, -3),
                (86, -5),
                (66, -11),
                (90, -12),
                (88, -8),
                (42, -15),
            ])
        ),
        (
            81,
            HashMap::from([
                (65, -3),
                (82, -1),
                (78, -3),
                (68, -2),
                (67, -12),
                (81, 8),
                (69, 2),
                (71, -6),
                (72, 1),
                (73, -7),
                (76, -4),
                (75, -2),
                (77, -3),
                (70, -11),
                (80, -2),
                (83, -4),
                (84, -5),
                (87, -11),
                (89, -10),
                (86, -6),
                (66, -2),
                (90, 6),
                (88, -4),
                (42, -15),
            ])
        ),
        (
            69,
            HashMap::from([
                (65, -2),
                (82, -8),
                (78, -1),
                (68, 3),
                (67, -12),
                (81, 2),
                (69, 7),
                (71, -3),
                (72, -4),
                (73, -5),
                (76, -8),
                (75, -4),
                (77, -6),
                (70, -12),
                (80, -5),
                (83, -4),
                (84, -5),
                (87, -15),
                (89, -8),
                (86, -6),
                (66, 2),
                (90, 6),
                (88, -4),
                (42, -15),
            ])
        ),
        (
            71,
            HashMap::from([
                (65, -1),
                (82, -8),
                (78, -2),
                (68, -3),
                (67, -8),
                (81, -6),
                (69, -3),
                (71, 6),
                (72, -8),
                (73, -9),
                (76, -9),
                (75, -6),
                (77, -7),
                (70, -8),
                (80, -5),
                (83, -1),
                (84, -5),
                (87, -13),
                (89, -12),
                (86, -5),
                (66, -2),
                (90, -4),
                (88, -4),
                (42, -15),
            ])
        ),
        (
            72,
            HashMap::from([
                (65, -6),
                (82, -1),
                (78, 1),
                (68, -3),
                (67, -7),
                (81, 1),
                (69, -4),
                (71, -8),
                (72, 9),
                (73, -8),
                (76, -5),
                (75, -5),
                (77, -9),
                (70, -5),
                (80, -3),
                (83, -5),
                (84, -6),
                (87, -6),
                (89, -3),
                (86, -6),
                (66, -1),
                (90, 0),
                (88, -4),
                (42, -15),
            ])
        ),
        (
            73,
            HashMap::from([
                (65, -4),
                (82, -5),
                (78, -4),
                (68, -6),
                (67, -5),
                (81, -7),
                (69, -5),
                (71, -9),
                (72, -8),
                (73, 8),
                (76, -1),
                (75, -5),
                (77, 0),
                (70, -2),
                (80, -7),
                (83, -6),
                (84, -2),
                (87, -12),
                (89, -5),
                (86, 2),
                (66, -5),
                (90, -5),
                (88, -4),
                (42, -15),
            ])
        ),
        (
            76,
            HashMap::from([
                (65, -5),
                (82, -8),
                (78, -6),
                (68, -11),
                (67, -13),
                (81, -4),
                (69, -8),
                (71, -9),
                (72, -5),
                (73, -1),
                (76, 7),
                (75, -7),
                (77, 1),
                (70, -2),
                (80, -6),
                (83, -7),
                (84, -6),
                (87, -5),
                (89, -6),
                (86, -2),
                (66, -8),
                (90, -6),
                (88, -5),
                (42, -15),
            ])
        ),
        (
            75,
            HashMap::from([
                (65, -6),
                (82, 1),
                (78, 0),
                (68, -4),
                (67, -12),
                (81, -2),
                (69, -4),
                (71, -6),
                (72, -5),
                (73, -5),
                (76, -7),
                (75, 6),
                (77, -1),
                (70, -12),
                (80, -6),
                (83, -3),
                (84, -2),
                (87, -10),
                (89, -8),
                (86, -8),
                (66, -2),
                (90, -3),
                (88, -4),
                (42, -15),
            ])
        ),
        (
            77,
            HashMap::from([
                (65, -4),
                (82, -3),
                (78, -7),
                (68, -9),
                (67, -12),
                (81, -3),
                (69, -6),
                (71, -7),
                (72, -9),
                (73, 0),
                (76, 1),
                (75, -1),
                (77, 11),
                (70, -3),
                (80, -7),
                (83, -5),
                (84, -3),
                (87, -11),
                (89, -10),
                (86, -1),
                (66, -8),
                (90, -4),
                (88, -4),
                (42, -15),
            ])
        ),
        (
            70,
            HashMap::from([
                (65, -7),
                (82, -8),
                (78, -8),
                (68, -13),
                (67, -11),
                (81, -11),
                (69, -12),
                (71, -8),
                (72, -5),
                (73, -2),
                (76, -2),
                (75, -12),
                (77, -3),
                (70, 9),
                (80, -9),
                (83, -6),
                (84, -8),
                (87, -4),
                (89, 2),
                (86, -7),
                (66, -9),
                (90, -12),
                (88, -7),
                (42, -15),
            ])
        ),
        (
            80,
            HashMap::from([
                (65, -1),
                (82, -3),
                (78, -5),
                (68, -7),
                (67, -7),
                (81, -2),
                (69, -5),
                (71, -5),
                (72, -3),
                (73, -7),
                (76, -6),
                (75, -6),
                (77, -7),
                (70, -9),
                (80, 8),
                (83, -1),
                (84, -3),
                (87, -12),
                (89, -12),
                (86, -5),
                (66, -6),
                (90, -3),
                (88, -4),
                (42, -15),
            ])
        ),
        (
            83,
            HashMap::from([
                (65, 0),
                (82, -2),
                (78, 0),
                (68, -3),
                (67, -2),
                (81, -4),
                (69, -4),
                (71, -1),
                (72, -5),
                (73, -6),
                (76, -7),
                (75, -3),
                (77, -5),
                (70, -6),
                (80, -1),
                (83, 6),
                (84, 1),
                (87, -4),
                (89, -6),
                (86, -5),
                (66, -1),
                (90, -4),
                (88, -2),
                (42, -15),
            ])
        ),
        (
            84,
            HashMap::from([
                (65, 0),
                (82, -5),
                (78, -1),
                (68, -4),
                (67, -7),
                (81, -5),
                (69, -5),
                (71, -5),
                (72, -6),
                (73, -2),
                (76, -6),
                (75, -2),
                (77, -3),
                (70, -8),
                (80, -3),
                (83, 1),
                (84, 7),
                (87, -11),
                (89, -6),
                (86, -2),
                (66, -2),
                (90, -5),
                (88, -3),
                (42, -15),
            ])
        ),
        (
            87,
            HashMap::from([
                (65, -12),
                (82, -1),
                (78, -7),
                (68, -13),
                (67, -14),
                (81, -11),
                (69, -15),
                (71, -13),
                (72, -6),
                (73, -12),
                (76, -5),
                (75, -10),
                (77, -11),
                (70, -4),
                (80, -12),
                (83, -4),
                (84, -11),
                (87, 13),
                (89, -4),
                (86, -14),
                (66, -9),
                (90, -13),
                (88, -9),
                (42, -15),
            ])
        ),
        (
            89,
            HashMap::from([
                (65, -7),
                (82, -9),
                (78, -4),
                (68, -10),
                (67, -3),
                (81, -10),
                (69, -8),
                (71, -12),
                (72, -3),
                (73, -5),
                (76, -6),
                (75, -8),
                (77, -10),
                (70, 2),
                (80, -12),
                (83, -6),
                (84, -6),
                (87, -4),
                (89, 10),
                (86, -6),
                (66, -6),
                (90, -8),
                (88, -7),
                (42, -15),
            ])
        ),
        (
            86,
            HashMap::from([
                (65, -2),
                (82, -7),
                (78, -7),
                (68, -7),
                (67, -5),
                (81, -6),
                (69, -6),
                (71, -5),
                (72, -6),
                (73, 2),
                (76, -2),
                (75, -8),
                (77, -1),
                (70, -7),
                (80, -5),
                (83, -5),
                (84, -2),
                (87, -14),
                (89, -6),
                (86, 7),
                (66, -7),
                (90, -6),
                (88, -4),
                (42, -15),
            ])
        ),
        (
            66,
            HashMap::from([
                (65, -3),
                (82, -6),
                (78, 6),
                (68, 6),
                (67, -11),
                (81, -2),
                (69, 2),
                (71, -2),
                (72, -1),
                (73, -5),
                (76, -8),
                (75, -2),
                (77, -8),
                (70, -9),
                (80, -6),
                (83, -1),
                (84, -2),
                (87, -9),
                (89, -6),
                (86, -7),
                (66, 6),
                (90, 1),
                (88, -4),
                (42, -15),
            ])
        ),
        (
            90,
            HashMap::from([
                (65, -2),
                (82, -3),
                (78, -2),
                (68, 2),
                (67, -12),
                (81, 6),
                (69, 6),
                (71, -4),
                (72, 0),
                (73, -5),
                (76, -6),
                (75, -3),
                (77, -4),
                (70, -12),
                (80, -3),
                (83, -4),
                (84, -5),
                (87, -13),
                (89, -8),
                (86, -6),
                (66, 1),
                (90, 6),
                (88, -4),
                (42, -15),
            ])
        ),
        (
            88,
            HashMap::from([
                (65, -3),
                (82, -5),
                (78, -3),
                (68, -5),
                (67, -8),
                (81, -4),
                (69, -4),
                (71, -4),
                (72, -4),
                (73, -4),
                (76, -5),
                (75, -4),
                (77, -4),
                (70, -7),
                (80, -4),
                (83, -2),
                (84, -3),
                (87, -9),
                (89, -7),
                (86, -4),
                (66, -4),
                (90, -4),
                (88, -4),
                (42, -15),
            ])
        ),
        (
            42,
            HashMap::from([
                (65, -15),
                (82, -15),
                (78, -15),
                (68, -15),
                (67, -15),
                (81, -15),
                (69, -15),
                (71, -15),
                (72, -15),
                (73, -15),
                (76, -15),
                (75, -15),
                (77, -15),
                (70, -15),
                (80, -15),
                (83, -15),
                (84, -15),
                (87, -15),
                (89, -15),
                (86, -15),
                (66, -15),
                (90, -15),
                (88, -15),
                (42, 1),
            ])
        ),
    ]);
}