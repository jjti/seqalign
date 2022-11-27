use super::Matrix;
use std::collections::HashMap;

lazy_static! {
    pub static ref MATRIX: Matrix = HashMap::from([
        (
            65,
            HashMap::from([
                (65, 3),
                (82, -3),
                (78, 0),
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
                (70, -5),
                (80, 1),
                (83, 2),
                (84, 2),
                (87, -8),
                (89, -5),
                (86, 0),
                (66, 0),
                (90, 0),
                (88, -1),
                (42, -10),
            ])
        ),
        (
            82,
            HashMap::from([
                (65, -3),
                (82, 8),
                (78, -1),
                (68, -3),
                (67, -5),
                (81, 1),
                (69, -2),
                (71, -4),
                (72, 2),
                (73, -3),
                (76, -4),
                (75, 4),
                (77, -1),
                (70, -6),
                (80, -1),
                (83, -1),
                (84, -2),
                (87, 2),
                (89, -6),
                (86, -4),
                (66, -2),
                (90, 0),
                (88, -2),
                (42, -10),
            ])
        ),
        (
            78,
            HashMap::from([
                (65, 0),
                (82, -1),
                (78, 4),
                (68, 3),
                (67, -5),
                (81, 0),
                (69, 2),
                (71, 0),
                (72, 2),
                (73, -3),
                (76, -4),
                (75, 1),
                (77, -3),
                (70, -5),
                (80, -1),
                (83, 1),
                (84, 0),
                (87, -5),
                (89, -3),
                (86, -3),
                (66, 3),
                (90, 1),
                (88, -1),
                (42, -10),
            ])
        ),
        (
            68,
            HashMap::from([
                (65, 0),
                (82, -3),
                (78, 3),
                (68, 6),
                (67, -7),
                (81, 2),
                (69, 5),
                (71, 0),
                (72, 0),
                (73, -4),
                (76, -6),
                (75, -1),
                (77, -4),
                (70, -8),
                (80, -2),
                (83, 0),
                (84, -1),
                (87, -9),
                (89, -6),
                (86, -4),
                (66, 5),
                (90, 4),
                (88, -1),
                (42, -10),
            ])
        ),
        (
            67,
            HashMap::from([
                (65, -3),
                (82, -5),
                (78, -5),
                (68, -7),
                (67, 13),
                (81, -8),
                (69, -8),
                (71, -5),
                (72, -5),
                (73, -3),
                (76, -9),
                (75, -8),
                (77, -7),
                (70, -6),
                (80, -4),
                (83, 0),
                (84, -3),
                (87, -10),
                (89, 0),
                (86, -3),
                (66, -6),
                (90, -8),
                (88, -4),
                (42, -10),
            ])
        ),
        (
            81,
            HashMap::from([
                (65, -1),
                (82, 1),
                (78, 0),
                (68, 2),
                (67, -8),
                (81, 6),
                (69, 3),
                (71, -2),
                (72, 4),
                (73, -3),
                (76, -2),
                (75, 0),
                (77, -1),
                (70, -7),
                (80, 0),
                (83, -1),
                (84, -2),
                (87, -7),
                (89, -6),
                (86, -3),
                (66, 1),
                (90, 5),
                (88, -1),
                (42, -10),
            ])
        ),
        (
            69,
            HashMap::from([
                (65, 0),
                (82, -2),
                (78, 2),
                (68, 5),
                (67, -8),
                (81, 3),
                (69, 6),
                (71, 0),
                (72, 0),
                (73, -3),
                (76, -5),
                (75, -1),
                (77, -3),
                (70, -8),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -10),
                (89, -6),
                (86, -3),
                (66, 3),
                (90, 5),
                (88, -1),
                (42, -10),
            ])
        ),
        (
            71,
            HashMap::from([
                (65, 1),
                (82, -4),
                (78, 0),
                (68, 0),
                (67, -5),
                (81, -2),
                (69, 0),
                (71, 6),
                (72, -4),
                (73, -4),
                (76, -6),
                (75, -3),
                (77, -4),
                (70, -6),
                (80, -2),
                (83, 1),
                (84, -1),
                (87, -9),
                (89, -7),
                (86, -2),
                (66, 0),
                (90, -1),
                (88, -2),
                (42, -10),
            ])
        ),
        (
            72,
            HashMap::from([
                (65, -3),
                (82, 2),
                (78, 2),
                (68, 0),
                (67, -5),
                (81, 4),
                (69, 0),
                (71, -4),
                (72, 9),
                (73, -4),
                (76, -3),
                (75, -1),
                (77, -4),
                (70, -3),
                (80, -1),
                (83, -2),
                (84, -2),
                (87, -4),
                (89, 0),
                (86, -3),
                (66, 1),
                (90, 2),
                (88, -1),
                (42, -10),
            ])
        ),
        (
            73,
            HashMap::from([
                (65, -1),
                (82, -3),
                (78, -3),
                (68, -4),
                (67, -3),
                (81, -3),
                (69, -3),
                (71, -4),
                (72, -4),
                (73, 7),
                (76, 2),
                (75, -3),
                (77, 2),
                (70, 1),
                (80, -3),
                (83, -2),
                (84, 0),
                (87, -7),
                (89, -2),
                (86, 5),
                (66, -3),
                (90, -3),
                (88, -1),
                (42, -10),
            ])
        ),
        (
            76,
            HashMap::from([
                (65, -3),
                (82, -4),
                (78, -4),
                (68, -6),
                (67, -9),
                (81, -2),
                (69, -5),
                (71, -6),
                (72, -3),
                (73, 2),
                (76, 7),
                (75, -4),
                (77, 4),
                (70, 1),
                (80, -4),
                (83, -4),
                (84, -3),
                (87, -3),
                (89, -2),
                (86, 2),
                (66, -5),
                (90, -4),
                (88, -2),
                (42, -10),
            ])
        ),
        (
            75,
            HashMap::from([
                (65, -2),
                (82, 4),
                (78, 1),
                (68, -1),
                (67, -8),
                (81, 0),
                (69, -1),
                (71, -3),
                (72, -1),
                (73, -3),
                (76, -4),
                (75, 6),
                (77, 1),
                (70, -8),
                (80, -2),
                (83, -1),
                (84, 0),
                (87, -5),
                (89, -6),
                (86, -4),
                (66, 0),
                (90, 0),
                (88, -2),
                (42, -10),
            ])
        ),
        (
            77,
            HashMap::from([
                (65, -2),
                (82, -1),
                (78, -3),
                (68, -4),
                (67, -7),
                (81, -1),
                (69, -3),
                (71, -4),
                (72, -4),
                (73, 2),
                (76, 4),
                (75, 1),
                (77, 10),
                (70, 0),
                (80, -3),
                (83, -2),
                (84, -1),
                (87, -6),
                (89, -4),
                (86, 2),
                (66, -4),
                (90, -2),
                (88, -1),
                (42, -10),
            ])
        ),
        (
            70,
            HashMap::from([
                (65, -5),
                (82, -6),
                (78, -5),
                (68, -8),
                (67, -6),
                (81, -7),
                (69, -8),
                (71, -6),
                (72, -3),
                (73, 1),
                (76, 1),
                (75, -8),
                (77, 0),
                (70, 10),
                (80, -6),
                (83, -4),
                (84, -5),
                (87, -1),
                (89, 7),
                (86, -2),
                (66, -6),
                (90, -7),
                (88, -4),
                (42, -10),
            ])
        ),
        (
            80,
            HashMap::from([
                (65, 1),
                (82, -1),
                (78, -1),
                (68, -2),
                (67, -4),
                (81, 0),
                (69, -1),
                (71, -2),
                (72, -1),
                (73, -3),
                (76, -4),
                (75, -2),
                (77, -3),
                (70, -6),
                (80, 8),
                (83, 1),
                (84, 0),
                (87, -8),
                (89, -7),
                (86, -2),
                (66, -2),
                (90, -1),
                (88, -1),
                (42, -10),
            ])
        ),
        (
            83,
            HashMap::from([
                (65, 2),
                (82, -1),
                (78, 1),
                (68, 0),
                (67, 0),
                (81, -1),
                (69, -1),
                (71, 1),
                (72, -2),
                (73, -2),
                (76, -4),
                (75, -1),
                (77, -2),
                (70, -4),
                (80, 1),
                (83, 3),
                (84, 2),
                (87, -3),
                (89, -4),
                (86, -2),
                (66, 1),
                (90, -1),
                (88, 0),
                (42, -10),
            ])
        ),
        (
            84,
            HashMap::from([
                (65, 2),
                (82, -2),
                (78, 0),
                (68, -1),
                (67, -3),
                (81, -2),
                (69, -1),
                (71, -1),
                (72, -2),
                (73, 0),
                (76, -3),
                (75, 0),
                (77, -1),
                (70, -5),
                (80, 0),
                (83, 2),
                (84, 5),
                (87, -7),
                (89, -4),
                (86, 0),
                (66, 0),
                (90, -1),
                (88, -1),
                (42, -10),
            ])
        ),
        (
            87,
            HashMap::from([
                (65, -8),
                (82, 2),
                (78, -5),
                (68, -9),
                (67, -10),
                (81, -7),
                (69, -10),
                (71, -9),
                (72, -4),
                (73, -7),
                (76, -3),
                (75, -5),
                (77, -6),
                (70, -1),
                (80, -8),
                (83, -3),
                (84, -7),
                (87, 18),
                (89, -1),
                (86, -9),
                (66, -7),
                (90, -8),
                (88, -6),
                (42, -10),
            ])
        ),
        (
            89,
            HashMap::from([
                (65, -5),
                (82, -6),
                (78, -3),
                (68, -6),
                (67, 0),
                (81, -6),
                (69, -6),
                (71, -7),
                (72, 0),
                (73, -2),
                (76, -2),
                (75, -6),
                (77, -4),
                (70, 7),
                (80, -7),
                (83, -4),
                (84, -4),
                (87, -1),
                (89, 12),
                (86, -4),
                (66, -4),
                (90, -6),
                (88, -4),
                (42, -10),
            ])
        ),
        (
            86,
            HashMap::from([
                (65, 0),
                (82, -4),
                (78, -3),
                (68, -4),
                (67, -3),
                (81, -3),
                (69, -3),
                (71, -2),
                (72, -3),
                (73, 5),
                (76, 2),
                (75, -4),
                (77, 2),
                (70, -2),
                (80, -2),
                (83, -2),
                (84, 0),
                (87, -9),
                (89, -4),
                (86, 6),
                (66, -3),
                (90, -3),
                (88, -1),
                (42, -10),
            ])
        ),
        (
            66,
            HashMap::from([
                (65, 0),
                (82, -2),
                (78, 3),
                (68, 5),
                (67, -6),
                (81, 1),
                (69, 3),
                (71, 0),
                (72, 1),
                (73, -3),
                (76, -5),
                (75, 0),
                (77, -4),
                (70, -6),
                (80, -2),
                (83, 1),
                (84, 0),
                (87, -7),
                (89, -4),
                (86, -3),
                (66, 4),
                (90, 3),
                (88, -1),
                (42, -10),
            ])
        ),
        (
            90,
            HashMap::from([
                (65, 0),
                (82, 0),
                (78, 1),
                (68, 4),
                (67, -8),
                (81, 5),
                (69, 5),
                (71, -1),
                (72, 2),
                (73, -3),
                (76, -4),
                (75, 0),
                (77, -2),
                (70, -7),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -8),
                (89, -6),
                (86, -3),
                (66, 3),
                (90, 5),
                (88, -1),
                (42, -10),
            ])
        ),
        (
            88,
            HashMap::from([
                (65, -1),
                (82, -2),
                (78, -1),
                (68, -1),
                (67, -4),
                (81, -1),
                (69, -1),
                (71, -2),
                (72, -1),
                (73, -1),
                (76, -2),
                (75, -2),
                (77, -1),
                (70, -4),
                (80, -1),
                (83, 0),
                (84, -1),
                (87, -6),
                (89, -4),
                (86, -1),
                (66, -1),
                (90, -1),
                (88, -2),
                (42, -10),
            ])
        ),
        (
            42,
            HashMap::from([
                (65, -10),
                (82, -10),
                (78, -10),
                (68, -10),
                (67, -10),
                (81, -10),
                (69, -10),
                (71, -10),
                (72, -10),
                (73, -10),
                (76, -10),
                (75, -10),
                (77, -10),
                (70, -10),
                (80, -10),
                (83, -10),
                (84, -10),
                (87, -10),
                (89, -10),
                (86, -10),
                (66, -10),
                (90, -10),
                (88, -10),
                (42, 1),
            ])
        ),
    ]);
}