use super::Matrix;
use std::collections::HashMap;

lazy_static! {
    pub static ref MATRIX: Matrix = HashMap::from([
        (
            65,
            HashMap::from([
                (65, 5),
                (82, -5),
                (78, -2),
                (68, -2),
                (67, -5),
                (81, -3),
                (69, -1),
                (71, -1),
                (72, -5),
                (73, -3),
                (76, -5),
                (75, -5),
                (77, -4),
                (70, -7),
                (80, 0),
                (83, 0),
                (84, 0),
                (87, -11),
                (89, -6),
                (86, -1),
                (66, -2),
                (90, -2),
                (88, -2),
                (42, -13),
            ])
        ),
        (
            82,
            HashMap::from([
                (65, -5),
                (82, 8),
                (78, -4),
                (68, -7),
                (67, -6),
                (81, 0),
                (69, -7),
                (71, -7),
                (72, 0),
                (73, -4),
                (76, -7),
                (75, 1),
                (77, -3),
                (70, -8),
                (80, -3),
                (83, -2),
                (84, -5),
                (87, -1),
                (89, -8),
                (86, -6),
                (66, -5),
                (90, -2),
                (88, -4),
                (42, -13),
            ])
        ),
        (
            78,
            HashMap::from([
                (65, -2),
                (82, -4),
                (78, 7),
                (68, 2),
                (67, -8),
                (81, -2),
                (69, -1),
                (71, -2),
                (72, 1),
                (73, -4),
                (76, -6),
                (75, 0),
                (77, -6),
                (70, -7),
                (80, -4),
                (83, 1),
                (84, -1),
                (87, -7),
                (89, -3),
                (86, -6),
                (66, 5),
                (90, -1),
                (88, -2),
                (42, -13),
            ])
        ),
        (
            68,
            HashMap::from([
                (65, -2),
                (82, -7),
                (78, 2),
                (68, 7),
                (67, -11),
                (81, -1),
                (69, 3),
                (71, -2),
                (72, -2),
                (73, -6),
                (76, -10),
                (75, -3),
                (77, -8),
                (70, -12),
                (80, -6),
                (83, -2),
                (84, -3),
                (87, -12),
                (89, -9),
                (86, -6),
                (66, 6),
                (90, 2),
                (88, -4),
                (42, -13),
            ])
        ),
        (
            67,
            HashMap::from([
                (65, -5),
                (82, -6),
                (78, -8),
                (68, -11),
                (67, 9),
                (81, -11),
                (69, -11),
                (71, -7),
                (72, -6),
                (73, -5),
                (76, -12),
                (75, -11),
                (77, -11),
                (70, -10),
                (80, -6),
                (83, -2),
                (84, -6),
                (87, -13),
                (89, -3),
                (86, -5),
                (66, -9),
                (90, -11),
                (88, -7),
                (42, -13),
            ])
        ),
        (
            81,
            HashMap::from([
                (65, -3),
                (82, 0),
                (78, -2),
                (68, -1),
                (67, -11),
                (81, 8),
                (69, 2),
                (71, -5),
                (72, 2),
                (73, -6),
                (76, -4),
                (75, -2),
                (77, -3),
                (70, -10),
                (80, -2),
                (83, -4),
                (84, -4),
                (87, -10),
                (89, -9),
                (86, -5),
                (66, -2),
                (90, 6),
                (88, -3),
                (42, -13),
            ])
        ),
        (
            69,
            HashMap::from([
                (65, -1),
                (82, -7),
                (78, -1),
                (68, 3),
                (67, -11),
                (81, 2),
                (69, 7),
                (71, -3),
                (72, -3),
                (73, -4),
                (76, -7),
                (75, -3),
                (77, -5),
                (70, -11),
                (80, -4),
                (83, -3),
                (84, -4),
                (87, -13),
                (89, -7),
                (86, -5),
                (66, 2),
                (90, 6),
                (88, -3),
                (42, -13),
            ])
        ),
        (
            71,
            HashMap::from([
                (65, -1),
                (82, -7),
                (78, -2),
                (68, -2),
                (67, -7),
                (81, -5),
                (69, -3),
                (71, 6),
                (72, -7),
                (73, -8),
                (76, -9),
                (75, -6),
                (77, -7),
                (70, -8),
                (80, -4),
                (83, -1),
                (84, -4),
                (87, -12),
                (89, -11),
                (86, -4),
                (66, -2),
                (90, -4),
                (88, -4),
                (42, -13),
            ])
        ),
        (
            72,
            HashMap::from([
                (65, -5),
                (82, 0),
                (78, 1),
                (68, -2),
                (67, -6),
                (81, 2),
                (69, -3),
                (71, -7),
                (72, 9),
                (73, -7),
                (76, -5),
                (75, -4),
                (77, -8),
                (70, -5),
                (80, -3),
                (83, -4),
                (84, -5),
                (87, -6),
                (89, -2),
                (86, -5),
                (66, 0),
                (90, 0),
                (88, -4),
                (42, -13),
            ])
        ),
        (
            73,
            HashMap::from([
                (65, -3),
                (82, -4),
                (78, -4),
                (68, -6),
                (67, -5),
                (81, -6),
                (69, -4),
                (71, -8),
                (72, -7),
                (73, 8),
                (76, 0),
                (75, -5),
                (77, 0),
                (70, -1),
                (80, -7),
                (83, -5),
                (84, -1),
                (87, -11),
                (89, -5),
                (86, 3),
                (66, -5),
                (90, -5),
                (88, -3),
                (42, -13),
            ])
        ),
        (
            76,
            HashMap::from([
                (65, -5),
                (82, -7),
                (78, -6),
                (68, -10),
                (67, -12),
                (81, -4),
                (69, -7),
                (71, -9),
                (72, -5),
                (73, 0),
                (76, 6),
                (75, -6),
                (77, 2),
                (70, -1),
                (80, -6),
                (83, -7),
                (84, -5),
                (87, -5),
                (89, -5),
                (86, -1),
                (66, -7),
                (90, -5),
                (88, -5),
                (42, -13),
            ])
        ),
        (
            75,
            HashMap::from([
                (65, -5),
                (82, 1),
                (78, 0),
                (68, -3),
                (67, -11),
                (81, -2),
                (69, -3),
                (71, -6),
                (72, -4),
                (73, -5),
                (76, -6),
                (75, 6),
                (77, -1),
                (70, -11),
                (80, -5),
                (83, -3),
                (84, -2),
                (87, -9),
                (89, -8),
                (86, -7),
                (66, -1),
                (90, -2),
                (88, -4),
                (42, -13),
            ])
        ),
        (
            77,
            HashMap::from([
                (65, -4),
                (82, -3),
                (78, -6),
                (68, -8),
                (67, -11),
                (81, -3),
                (69, -5),
                (71, -7),
                (72, -8),
                (73, 0),
                (76, 2),
                (75, -1),
                (77, 10),
                (70, -3),
                (80, -6),
                (83, -4),
                (84, -3),
                (87, -10),
                (89, -8),
                (86, 0),
                (66, -7),
                (90, -4),
                (88, -4),
                (42, -13),
            ])
        ),
        (
            70,
            HashMap::from([
                (65, -7),
                (82, -8),
                (78, -7),
                (68, -12),
                (67, -10),
                (81, -10),
                (69, -11),
                (71, -8),
                (72, -5),
                (73, -1),
                (76, -1),
                (75, -11),
                (77, -3),
                (70, 9),
                (80, -8),
                (83, -5),
                (84, -7),
                (87, -3),
                (89, 3),
                (86, -6),
                (66, -9),
                (90, -11),
                (88, -6),
                (42, -13),
            ])
        ),
        (
            80,
            HashMap::from([
                (65, 0),
                (82, -3),
                (78, -4),
                (68, -6),
                (67, -6),
                (81, -2),
                (69, -4),
                (71, -4),
                (72, -3),
                (73, -7),
                (76, -6),
                (75, -5),
                (77, -6),
                (70, -8),
                (80, 8),
                (83, -1),
                (84, -3),
                (87, -11),
                (89, -11),
                (86, -4),
                (66, -5),
                (90, -3),
                (88, -4),
                (42, -13),
            ])
        ),
        (
            83,
            HashMap::from([
                (65, 0),
                (82, -2),
                (78, 1),
                (68, -2),
                (67, -2),
                (81, -4),
                (69, -3),
                (71, -1),
                (72, -4),
                (73, -5),
                (76, -7),
                (75, -3),
                (77, -4),
                (70, -5),
                (80, -1),
                (83, 6),
                (84, 1),
                (87, -4),
                (89, -5),
                (86, -4),
                (66, -1),
                (90, -3),
                (88, -2),
                (42, -13),
            ])
        ),
        (
            84,
            HashMap::from([
                (65, 0),
                (82, -5),
                (78, -1),
                (68, -3),
                (67, -6),
                (81, -4),
                (69, -4),
                (71, -4),
                (72, -5),
                (73, -1),
                (76, -5),
                (75, -2),
                (77, -3),
                (70, -7),
                (80, -3),
                (83, 1),
                (84, 6),
                (87, -10),
                (89, -5),
                (86, -2),
                (66, -2),
                (90, -4),
                (88, -2),
                (42, -13),
            ])
        ),
        (
            87,
            HashMap::from([
                (65, -11),
                (82, -1),
                (78, -7),
                (68, -12),
                (67, -13),
                (81, -10),
                (69, -13),
                (71, -12),
                (72, -6),
                (73, -11),
                (76, -5),
                (75, -9),
                (77, -10),
                (70, -3),
                (80, -11),
                (83, -4),
                (84, -10),
                (87, 13),
                (89, -4),
                (86, -12),
                (66, -8),
                (90, -11),
                (88, -9),
                (42, -13),
            ])
        ),
        (
            89,
            HashMap::from([
                (65, -6),
                (82, -8),
                (78, -3),
                (68, -9),
                (67, -3),
                (81, -9),
                (69, -7),
                (71, -11),
                (72, -2),
                (73, -5),
                (76, -5),
                (75, -8),
                (77, -8),
                (70, 3),
                (80, -11),
                (83, -5),
                (84, -5),
                (87, -4),
                (89, 9),
                (86, -6),
                (66, -5),
                (90, -8),
                (88, -6),
                (42, -13),
            ])
        ),
        (
            86,
            HashMap::from([
                (65, -1),
                (82, -6),
                (78, -6),
                (68, -6),
                (67, -5),
                (81, -5),
                (69, -5),
                (71, -4),
                (72, -5),
                (73, 3),
                (76, -1),
                (75, -7),
                (77, 0),
                (70, -6),
                (80, -4),
                (83, -4),
                (84, -2),
                (87, -12),
                (89, -6),
                (86, 7),
                (66, -6),
                (90, -5),
                (88, -3),
                (42, -13),
            ])
        ),
        (
            66,
            HashMap::from([
                (65, -2),
                (82, -5),
                (78, 5),
                (68, 6),
                (67, -9),
                (81, -2),
                (69, 2),
                (71, -2),
                (72, 0),
                (73, -5),
                (76, -7),
                (75, -1),
                (77, -7),
                (70, -9),
                (80, -5),
                (83, -1),
                (84, -2),
                (87, -8),
                (89, -5),
                (86, -6),
                (66, 5),
                (90, 1),
                (88, -3),
                (42, -13),
            ])
        ),
        (
            90,
            HashMap::from([
                (65, -2),
                (82, -2),
                (78, -1),
                (68, 2),
                (67, -11),
                (81, 6),
                (69, 6),
                (71, -4),
                (72, 0),
                (73, -5),
                (76, -5),
                (75, -2),
                (77, -4),
                (70, -11),
                (80, -3),
                (83, -3),
                (84, -4),
                (87, -11),
                (89, -8),
                (86, -5),
                (66, 1),
                (90, 6),
                (88, -3),
                (42, -13),
            ])
        ),
        (
            88,
            HashMap::from([
                (65, -2),
                (82, -4),
                (78, -2),
                (68, -4),
                (67, -7),
                (81, -3),
                (69, -3),
                (71, -4),
                (72, -4),
                (73, -3),
                (76, -5),
                (75, -4),
                (77, -4),
                (70, -6),
                (80, -4),
                (83, -2),
                (84, -2),
                (87, -9),
                (89, -6),
                (86, -3),
                (66, -3),
                (90, -3),
                (88, -4),
                (42, -13),
            ])
        ),
        (
            42,
            HashMap::from([
                (65, -13),
                (82, -13),
                (78, -13),
                (68, -13),
                (67, -13),
                (81, -13),
                (69, -13),
                (71, -13),
                (72, -13),
                (73, -13),
                (76, -13),
                (75, -13),
                (77, -13),
                (70, -13),
                (80, -13),
                (83, -13),
                (84, -13),
                (87, -13),
                (89, -13),
                (86, -13),
                (66, -13),
                (90, -13),
                (88, -13),
                (42, 1),
            ])
        ),
    ]);
}
