use super::Matrix;
use std::collections::HashMap;

lazy_static! {
    pub static ref MATRIX: Matrix = HashMap::from([
        (
            65,
            HashMap::from([
                (65, 2),
                (82, -1),
                (78, 0),
                (68, 0),
                (67, -2),
                (81, 0),
                (69, 0),
                (71, 2),
                (72, -1),
                (73, 0),
                (76, -2),
                (75, -1),
                (77, -1),
                (70, -4),
                (80, 1),
                (83, 1),
                (84, 1),
                (87, -6),
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
                (65, -1),
                (82, 7),
                (78, 0),
                (68, -1),
                (67, -4),
                (81, 2),
                (69, -1),
                (71, -2),
                (72, 2),
                (73, -2),
                (76, -3),
                (75, 4),
                (77, 0),
                (70, -5),
                (80, 0),
                (83, 0),
                (84, -1),
                (87, 3),
                (89, -5),
                (86, -3),
                (66, 0),
                (90, 0),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            78,
            HashMap::from([
                (65, 0),
                (82, 0),
                (78, 2),
                (68, 2),
                (67, -4),
                (81, 1),
                (69, 2),
                (71, 1),
                (72, 2),
                (73, -2),
                (76, -3),
                (75, 1),
                (77, -2),
                (70, -4),
                (80, 0),
                (83, 1),
                (84, 0),
                (87, -5),
                (89, -2),
                (86, -2),
                (66, 2),
                (90, 1),
                (88, 0),
                (42, -9),
            ])
        ),
        (
            68,
            HashMap::from([
                (65, 0),
                (82, -1),
                (78, 2),
                (68, 4),
                (67, -6),
                (81, 2),
                (69, 4),
                (71, 1),
                (72, 1),
                (73, -2),
                (76, -4),
                (75, 0),
                (77, -3),
                (70, -6),
                (80, -1),
                (83, 0),
                (84, 0),
                (87, -7),
                (89, -5),
                (86, -2),
                (66, 3),
                (90, 3),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            67,
            HashMap::from([
                (65, -2),
                (82, -4),
                (78, -4),
                (68, -6),
                (67, 15),
                (81, -6),
                (69, -6),
                (71, -4),
                (72, -4),
                (73, -3),
                (76, -7),
                (75, -6),
                (77, -6),
                (70, -5),
                (80, -3),
                (83, 0),
                (84, -2),
                (87, -9),
                (89, 1),
                (86, -2),
                (66, -5),
                (90, -6),
                (88, -3),
                (42, -9),
            ])
        ),
        (
            81,
            HashMap::from([
                (65, 0),
                (82, 2),
                (78, 1),
                (68, 2),
                (67, -6),
                (81, 4),
                (69, 3),
                (71, -1),
                (72, 3),
                (73, -2),
                (76, -2),
                (75, 1),
                (77, -1),
                (70, -5),
                (80, 0),
                (83, 0),
                (84, -1),
                (87, -5),
                (89, -4),
                (86, -2),
                (66, 2),
                (90, 3),
                (88, 0),
                (42, -9),
            ])
        ),
        (
            69,
            HashMap::from([
                (65, 0),
                (82, -1),
                (78, 2),
                (68, 4),
                (67, -6),
                (81, 3),
                (69, 4),
                (71, 0),
                (72, 1),
                (73, -2),
                (76, -4),
                (75, 0),
                (77, -2),
                (70, -6),
                (80, 0),
                (83, 0),
                (84, 0),
                (87, -8),
                (89, -5),
                (86, -2),
                (66, 3),
                (90, 3),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            71,
            HashMap::from([
                (65, 2),
                (82, -2),
                (78, 1),
                (68, 1),
                (67, -4),
                (81, -1),
                (69, 0),
                (71, 5),
                (72, -2),
                (73, -3),
                (76, -4),
                (75, -2),
                (77, -3),
                (70, -5),
                (80, 0),
                (83, 1),
                (84, 0),
                (87, -8),
                (89, -6),
                (86, -1),
                (66, 1),
                (90, 0),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            72,
            HashMap::from([
                (65, -1),
                (82, 2),
                (78, 2),
                (68, 1),
                (67, -4),
                (81, 3),
                (69, 1),
                (71, -2),
                (72, 7),
                (73, -2),
                (76, -2),
                (75, 0),
                (77, -2),
                (70, -2),
                (80, 0),
                (83, -1),
                (84, -1),
                (87, -3),
                (89, 0),
                (86, -2),
                (66, 1),
                (90, 2),
                (88, 0),
                (42, -9),
            ])
        ),
        (
            73,
            HashMap::from([
                (65, 0),
                (82, -2),
                (78, -2),
                (68, -2),
                (67, -3),
                (81, -2),
                (69, -2),
                (71, -3),
                (72, -2),
                (73, 5),
                (76, 3),
                (75, -2),
                (77, 3),
                (70, 1),
                (80, -2),
                (83, -1),
                (84, 0),
                (87, -6),
                (89, -1),
                (86, 4),
                (66, -2),
                (90, -2),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            76,
            HashMap::from([
                (65, -2),
                (82, -3),
                (78, -3),
                (68, -4),
                (67, -7),
                (81, -2),
                (69, -4),
                (71, -4),
                (72, -2),
                (73, 3),
                (76, 7),
                (75, -3),
                (77, 4),
                (70, 3),
                (80, -3),
                (83, -3),
                (84, -2),
                (87, -2),
                (89, 0),
                (86, 2),
                (66, -4),
                (90, -3),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            75,
            HashMap::from([
                (65, -1),
                (82, 4),
                (78, 1),
                (68, 0),
                (67, -6),
                (81, 1),
                (69, 0),
                (71, -2),
                (72, 0),
                (73, -2),
                (76, -3),
                (75, 5),
                (77, 0),
                (70, -6),
                (80, -1),
                (83, 0),
                (84, 0),
                (87, -4),
                (89, -5),
                (86, -2),
                (66, 1),
                (90, 1),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            77,
            HashMap::from([
                (65, -1),
                (82, 0),
                (78, -2),
                (68, -3),
                (67, -6),
                (81, -1),
                (69, -2),
                (71, -3),
                (72, -2),
                (73, 3),
                (76, 4),
                (75, 0),
                (77, 6),
                (70, 1),
                (80, -2),
                (83, -2),
                (84, -1),
                (87, -5),
                (89, -2),
                (86, 2),
                (66, -2),
                (90, -2),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            70,
            HashMap::from([
                (65, -4),
                (82, -5),
                (78, -4),
                (68, -6),
                (67, -5),
                (81, -5),
                (69, -6),
                (71, -5),
                (72, -2),
                (73, 1),
                (76, 3),
                (75, -6),
                (77, 1),
                (70, 11),
                (80, -5),
                (83, -4),
                (84, -3),
                (87, 1),
                (89, 9),
                (86, -1),
                (66, -5),
                (90, -5),
                (88, -2),
                (42, -9),
            ])
        ),
        (
            80,
            HashMap::from([
                (65, 1),
                (82, 0),
                (78, 0),
                (68, -1),
                (67, -3),
                (81, 0),
                (69, 0),
                (71, 0),
                (72, 0),
                (73, -2),
                (76, -3),
                (75, -1),
                (77, -2),
                (70, -5),
                (80, 6),
                (83, 1),
                (84, 1),
                (87, -6),
                (89, -5),
                (86, -1),
                (66, 0),
                (90, 0),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            83,
            HashMap::from([
                (65, 1),
                (82, 0),
                (78, 1),
                (68, 0),
                (67, 0),
                (81, 0),
                (69, 0),
                (71, 1),
                (72, -1),
                (73, -1),
                (76, -3),
                (75, 0),
                (77, -2),
                (70, -4),
                (80, 1),
                (83, 1),
                (84, 1),
                (87, -3),
                (89, -3),
                (86, -1),
                (66, 1),
                (90, 0),
                (88, 0),
                (42, -9),
            ])
        ),
        (
            84,
            HashMap::from([
                (65, 1),
                (82, -1),
                (78, 0),
                (68, 0),
                (67, -2),
                (81, -1),
                (69, 0),
                (71, 0),
                (72, -1),
                (73, 0),
                (76, -2),
                (75, 0),
                (77, -1),
                (70, -3),
                (80, 1),
                (83, 1),
                (84, 2),
                (87, -6),
                (89, -3),
                (86, 0),
                (66, 0),
                (90, 0),
                (88, 0),
                (42, -9),
            ])
        ),
        (
            87,
            HashMap::from([
                (65, -6),
                (82, 3),
                (78, -5),
                (68, -7),
                (67, -9),
                (81, -5),
                (69, -8),
                (71, -8),
                (72, -3),
                (73, -6),
                (76, -2),
                (75, -4),
                (77, -5),
                (70, 1),
                (80, -6),
                (83, -3),
                (84, -6),
                (87, 22),
                (89, 0),
                (86, -7),
                (66, -6),
                (90, -6),
                (88, -4),
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
                (67, 1),
                (81, -4),
                (69, -5),
                (71, -6),
                (72, 0),
                (73, -1),
                (76, 0),
                (75, -5),
                (77, -2),
                (70, 9),
                (80, -5),
                (83, -3),
                (84, -3),
                (87, 0),
                (89, 12),
                (86, -3),
                (66, -4),
                (90, -5),
                (88, -2),
                (42, -9),
            ])
        ),
        (
            86,
            HashMap::from([
                (65, 0),
                (82, -3),
                (78, -2),
                (68, -2),
                (67, -2),
                (81, -2),
                (69, -2),
                (71, -1),
                (72, -2),
                (73, 4),
                (76, 2),
                (75, -2),
                (77, 2),
                (70, -1),
                (80, -1),
                (83, -1),
                (84, 0),
                (87, -7),
                (89, -3),
                (86, 5),
                (66, -2),
                (90, -2),
                (88, 0),
                (42, -9),
            ])
        ),
        (
            66,
            HashMap::from([
                (65, 0),
                (82, 0),
                (78, 2),
                (68, 3),
                (67, -5),
                (81, 2),
                (69, 3),
                (71, 1),
                (72, 1),
                (73, -2),
                (76, -4),
                (75, 1),
                (77, -2),
                (70, -5),
                (80, 0),
                (83, 1),
                (84, 0),
                (87, -6),
                (89, -4),
                (86, -2),
                (66, 3),
                (90, 2),
                (88, 0),
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
                (67, -6),
                (81, 3),
                (69, 3),
                (71, 0),
                (72, 2),
                (73, -2),
                (76, -3),
                (75, 1),
                (77, -2),
                (70, -5),
                (80, 0),
                (83, 0),
                (84, 0),
                (87, -6),
                (89, -5),
                (86, -2),
                (66, 2),
                (90, 3),
                (88, -1),
                (42, -9),
            ])
        ),
        (
            88,
            HashMap::from([
                (65, 0),
                (82, -1),
                (78, 0),
                (68, -1),
                (67, -3),
                (81, 0),
                (69, -1),
                (71, -1),
                (72, 0),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -2),
                (80, -1),
                (83, 0),
                (84, 0),
                (87, -4),
                (89, -2),
                (86, 0),
                (66, 0),
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
