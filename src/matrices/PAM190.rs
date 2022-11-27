use super::Matrix;
use std::collections::HashMap;

lazy_static! {
    pub static ref MATRIX: Matrix = HashMap::from([
        (
            65,
            HashMap::from([
                (65, 3),
                (82, -2),
                (78, 0),
                (68, 0),
                (67, -3),
                (81, -1),
                (69, 0),
                (71, 1),
                (72, -2),
                (73, -1),
                (76, -3),
                (75, -2),
                (77, -2),
                (70, -5),
                (80, 1),
                (83, 1),
                (84, 2),
                (87, -7),
                (89, -4),
                (86, 0),
                (66, 0),
                (90, 0),
                (88, 0),
                (42, -9),
            ])
        ),
        (
            82,
            HashMap::from([
                (65, -2),
                (82, 8),
                (78, -1),
                (68, -2),
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
                (80, 0),
                (83, -1),
                (84, -2),
                (87, 2),
                (89, -5),
                (86, -3),
                (66, -1),
                (90, 0),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            78,
            HashMap::from([
                (65, 0),
                (82, -1),
                (78, 3),
                (68, 3),
                (67, -5),
                (81, 1),
                (69, 2),
                (71, 0),
                (72, 2),
                (73, -2),
                (76, -4),
                (75, 1),
                (77, -3),
                (70, -4),
                (80, -1),
                (83, 1),
                (84, 0),
                (87, -5),
                (89, -2),
                (86, -3),
                (66, 3),
                (90, 1),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            68,
            HashMap::from([
                (65, 0),
                (82, -2),
                (78, 3),
                (68, 5),
                (67, -7),
                (81, 2),
                (69, 4),
                (71, 0),
                (72, 0),
                (73, -3),
                (76, -5),
                (75, 0),
                (77, -4),
                (70, -7),
                (80, -2),
                (83, 0),
                (84, -1),
                (87, -8),
                (89, -5),
                (86, -3),
                (66, 4),
                (90, 3),
                (88, -1),
                (42, -9),
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
                (81, -7),
                (69, -7),
                (71, -4),
                (72, -4),
                (73, -3),
                (76, -8),
                (75, -7),
                (77, -7),
                (70, -6),
                (80, -4),
                (83, 0),
                (84, -3),
                (87, -9),
                (89, 0),
                (86, -3),
                (66, -6),
                (90, -7),
                (88, -4),
                (42, -9),
            ])
        ),
        (
            81,
            HashMap::from([
                (65, -1),
                (82, 1),
                (78, 1),
                (68, 2),
                (67, -7),
                (81, 6),
                (69, 3),
                (71, -2),
                (72, 3),
                (73, -3),
                (76, -2),
                (75, 1),
                (77, -1),
                (70, -6),
                (80, 0),
                (83, -1),
                (84, -1),
                (87, -6),
                (89, -5),
                (86, -3),
                (66, 1),
                (90, 4),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            69,
            HashMap::from([
                (65, 0),
                (82, -2),
                (78, 2),
                (68, 4),
                (67, -7),
                (81, 3),
                (69, 5),
                (71, 0),
                (72, 0),
                (73, -3),
                (76, -4),
                (75, -1),
                (77, -3),
                (70, -7),
                (80, -1),
                (83, 0),
                (84, -1),
                (87, -9),
                (89, -5),
                (86, -3),
                (66, 3),
                (90, 4),
                (88, -1),
                (42, -9),
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
                (81, -2),
                (69, 0),
                (71, 6),
                (72, -3),
                (73, -4),
                (76, -5),
                (75, -3),
                (77, -4),
                (70, -6),
                (80, -1),
                (83, 1),
                (84, -1),
                (87, -9),
                (89, -7),
                (86, -2),
                (66, 0),
                (90, -1),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            72,
            HashMap::from([
                (65, -2),
                (82, 2),
                (78, 2),
                (68, 0),
                (67, -4),
                (81, 3),
                (69, 0),
                (71, -3),
                (72, 8),
                (73, -3),
                (76, -3),
                (75, -1),
                (77, -3),
                (70, -2),
                (80, -1),
                (83, -1),
                (84, -2),
                (87, -3),
                (89, 0),
                (86, -3),
                (66, 1),
                (90, 2),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            73,
            HashMap::from([
                (65, -1),
                (82, -3),
                (78, -2),
                (68, -3),
                (67, -3),
                (81, -3),
                (69, -3),
                (71, -4),
                (72, -3),
                (73, 6),
                (76, 2),
                (75, -3),
                (77, 2),
                (70, 1),
                (80, -3),
                (83, -2),
                (84, 0),
                (87, -7),
                (89, -2),
                (86, 4),
                (66, -3),
                (90, -3),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            76,
            HashMap::from([
                (65, -3),
                (82, -4),
                (78, -4),
                (68, -5),
                (67, -8),
                (81, -2),
                (69, -4),
                (71, -5),
                (72, -3),
                (73, 2),
                (76, 7),
                (75, -4),
                (77, 4),
                (70, 2),
                (80, -3),
                (83, -4),
                (84, -2),
                (87, -3),
                (89, -2),
                (86, 2),
                (66, -5),
                (90, -3),
                (88, -2),
                (42, -9),
            ])
        ),
        (
            75,
            HashMap::from([
                (65, -2),
                (82, 4),
                (78, 1),
                (68, 0),
                (67, -7),
                (81, 1),
                (69, -1),
                (71, -3),
                (72, -1),
                (73, -3),
                (76, -4),
                (75, 6),
                (77, 1),
                (70, -7),
                (80, -2),
                (83, 0),
                (84, 0),
                (87, -5),
                (89, -6),
                (86, -3),
                (66, 0),
                (90, 0),
                (88, -1),
                (42, -9),
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
                (72, -3),
                (73, 2),
                (76, 4),
                (75, 1),
                (77, 9),
                (70, 0),
                (80, -3),
                (83, -2),
                (84, -1),
                (87, -6),
                (89, -4),
                (86, 2),
                (66, -3),
                (90, -2),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            70,
            HashMap::from([
                (65, -5),
                (82, -6),
                (78, -4),
                (68, -7),
                (67, -6),
                (81, -6),
                (69, -7),
                (71, -6),
                (72, -2),
                (73, 1),
                (76, 2),
                (75, -7),
                (77, 0),
                (70, 10),
                (80, -6),
                (83, -4),
                (84, -4),
                (87, 0),
                (89, 7),
                (86, -2),
                (66, -6),
                (90, -7),
                (88, -3),
                (42, -9),
            ])
        ),
        (
            80,
            HashMap::from([
                (65, 1),
                (82, 0),
                (78, -1),
                (68, -2),
                (67, -4),
                (81, 0),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -3),
                (76, -3),
                (75, -2),
                (77, -3),
                (70, -6),
                (80, 7),
                (83, 1),
                (84, 0),
                (87, -7),
                (89, -6),
                (86, -2),
                (66, -1),
                (90, -1),
                (88, -1),
                (42, -9),
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
                (81, -1),
                (69, 0),
                (71, 1),
                (72, -1),
                (73, -2),
                (76, -4),
                (75, 0),
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
                (42, -9),
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
                (81, -1),
                (69, -1),
                (71, -1),
                (72, -2),
                (73, 0),
                (76, -2),
                (75, 0),
                (77, -1),
                (70, -4),
                (80, 0),
                (83, 2),
                (84, 4),
                (87, -6),
                (89, -3),
                (86, 0),
                (66, 0),
                (90, -1),
                (88, 0),
                (42, -9),
            ])
        ),
        (
            87,
            HashMap::from([
                (65, -7),
                (82, 2),
                (78, -5),
                (68, -8),
                (67, -9),
                (81, -6),
                (69, -9),
                (71, -9),
                (72, -3),
                (73, -7),
                (76, -3),
                (75, -5),
                (77, -6),
                (70, 0),
                (80, -7),
                (83, -3),
                (84, -6),
                (87, 18),
                (89, -1),
                (86, -8),
                (66, -6),
                (90, -7),
                (88, -5),
                (42, -9),
            ])
        ),
        (
            89,
            HashMap::from([
                (65, -4),
                (82, -5),
                (78, -2),
                (68, -5),
                (67, 0),
                (81, -5),
                (69, -5),
                (71, -7),
                (72, 0),
                (73, -2),
                (76, -2),
                (75, -6),
                (77, -4),
                (70, 7),
                (80, -6),
                (83, -4),
                (84, -3),
                (87, -1),
                (89, 11),
                (86, -3),
                (66, -4),
                (90, -5),
                (88, -3),
                (42, -9),
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
                (73, 4),
                (76, 2),
                (75, -3),
                (77, 2),
                (70, -2),
                (80, -2),
                (83, -2),
                (84, 0),
                (87, -8),
                (89, -3),
                (86, 6),
                (66, -3),
                (90, -3),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            66,
            HashMap::from([
                (65, 0),
                (82, -1),
                (78, 3),
                (68, 4),
                (67, -6),
                (81, 1),
                (69, 3),
                (71, 0),
                (72, 1),
                (73, -3),
                (76, -5),
                (75, 0),
                (77, -3),
                (70, -6),
                (80, -1),
                (83, 1),
                (84, 0),
                (87, -6),
                (89, -4),
                (86, -3),
                (66, 4),
                (90, 2),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            90,
            HashMap::from([
                (65, 0),
                (82, 0),
                (78, 1),
                (68, 3),
                (67, -7),
                (81, 4),
                (69, 4),
                (71, -1),
                (72, 2),
                (73, -3),
                (76, -3),
                (75, 0),
                (77, -2),
                (70, -7),
                (80, -1),
                (83, -1),
                (84, -1),
                (87, -7),
                (89, -5),
                (86, -3),
                (66, 2),
                (90, 4),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            88,
            HashMap::from([
                (65, 0),
                (82, -1),
                (78, -1),
                (68, -1),
                (67, -4),
                (81, -1),
                (69, -1),
                (71, -1),
                (72, -1),
                (73, -1),
                (76, -2),
                (75, -1),
                (77, -1),
                (70, -3),
                (80, -1),
                (83, 0),
                (84, 0),
                (87, -5),
                (89, -3),
                (86, -1),
                (66, -1),
                (90, -1),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            42,
            HashMap::from([
                (65, -9),
                (82, -9),
                (78, -9),
                (68, -9),
                (67, -9),
                (81, -9),
                (69, -9),
                (71, -9),
                (72, -9),
                (73, -9),
                (76, -9),
                (75, -9),
                (77, -9),
                (70, -9),
                (80, -9),
                (83, -9),
                (84, -9),
                (87, -9),
                (89, -9),
                (86, -9),
                (66, -9),
                (90, -9),
                (88, -9),
                (42, 1),
            ])
        ),
    ]);
}
