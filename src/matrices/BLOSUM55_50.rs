use super::Matrix;
use std::collections::HashMap;

lazy_static! {
    pub static ref MATRIX: Matrix = HashMap::from([
        (
            65,
            HashMap::from([
                (65, 5),
                (82, -2),
                (78, -2),
                (68, -2),
                (67, 0),
                (81, -1),
                (69, -1),
                (71, 0),
                (72, -2),
                (73, -2),
                (76, -2),
                (75, -1),
                (77, -1),
                (70, -3),
                (80, -1),
                (83, 2),
                (84, 0),
                (87, -4),
                (89, -2),
                (86, 0),
                (66, -2),
                (90, -1),
                (88, -1),
                (42, -5),
            ])
        ),
        (
            82,
            HashMap::from([
                (65, -2),
                (82, 8),
                (78, -1),
                (68, -2),
                (67, -4),
                (81, 1),
                (69, 0),
                (71, -3),
                (72, 0),
                (73, -4),
                (76, -3),
                (75, 3),
                (77, -2),
                (70, -3),
                (80, -3),
                (83, -1),
                (84, -1),
                (87, -3),
                (89, -2),
                (86, -3),
                (66, -1),
                (90, 0),
                (88, -1),
                (42, -5),
            ])
        ),
        (
            78,
            HashMap::from([
                (65, -2),
                (82, -1),
                (78, 8),
                (68, 2),
                (67, -3),
                (81, 0),
                (69, 0),
                (71, 0),
                (72, 1),
                (73, -4),
                (76, -4),
                (75, 0),
                (77, -3),
                (70, -4),
                (80, -2),
                (83, 1),
                (84, 0),
                (87, -5),
                (89, -2),
                (86, -4),
                (66, 4),
                (90, 0),
                (88, -1),
                (42, -5),
            ])
        ),
        (
            68,
            HashMap::from([
                (65, -2),
                (82, -2),
                (78, 2),
                (68, 8),
                (67, -4),
                (81, 0),
                (69, 2),
                (71, -2),
                (72, -1),
                (73, -4),
                (76, -5),
                (75, -1),
                (77, -4),
                (70, -5),
                (80, -2),
                (83, 0),
                (84, -1),
                (87, -5),
                (89, -3),
                (86, -4),
                (66, 5),
                (90, 1),
                (88, -2),
                (42, -5),
            ])
        ),
        (
            67,
            HashMap::from([
                (65, 0),
                (82, -4),
                (78, -3),
                (68, -4),
                (67, 13),
                (81, -4),
                (69, -4),
                (71, -3),
                (72, -4),
                (73, -2),
                (76, -2),
                (75, -4),
                (77, -2),
                (70, -3),
                (80, -3),
                (83, -1),
                (84, -1),
                (87, -4),
                (89, -3),
                (86, -1),
                (66, -4),
                (90, -4),
                (88, -2),
                (42, -5),
            ])
        ),
        (
            81,
            HashMap::from([
                (65, -1),
                (82, 1),
                (78, 0),
                (68, 0),
                (67, -4),
                (81, 7),
                (69, 2),
                (71, -2),
                (72, 1),
                (73, -4),
                (76, -3),
                (75, 2),
                (77, 0),
                (70, -4),
                (80, -1),
                (83, 0),
                (84, -1),
                (87, -2),
                (89, -1),
                (86, -3),
                (66, 0),
                (90, 4),
                (88, -1),
                (42, -5),
            ])
        ),
        (
            69,
            HashMap::from([
                (65, -1),
                (82, 0),
                (78, 0),
                (68, 2),
                (67, -4),
                (81, 2),
                (69, 7),
                (71, -3),
                (72, -1),
                (73, -4),
                (76, -4),
                (75, 1),
                (77, -3),
                (70, -4),
                (80, -1),
                (83, 0),
                (84, -1),
                (87, -3),
                (89, -2),
                (86, -3),
                (66, 1),
                (90, 5),
                (88, -1),
                (42, -5),
            ])
        ),
        (
            71,
            HashMap::from([
                (65, 0),
                (82, -3),
                (78, 0),
                (68, -2),
                (67, -3),
                (81, -2),
                (69, -3),
                (71, 8),
                (72, -2),
                (73, -5),
                (76, -5),
                (75, -2),
                (77, -3),
                (70, -4),
                (80, -3),
                (83, 0),
                (84, -2),
                (87, -3),
                (89, -4),
                (86, -4),
                (66, -1),
                (90, -3),
                (88, -2),
                (42, -5),
            ])
        ),
        (
            72,
            HashMap::from([
                (65, -2),
                (82, 0),
                (78, 1),
                (68, -1),
                (67, -4),
                (81, 1),
                (69, -1),
                (71, -2),
                (72, 11),
                (73, -4),
                (76, -3),
                (75, 0),
                (77, -2),
                (70, -1),
                (80, -3),
                (83, -1),
                (84, -2),
                (87, -3),
                (89, 2),
                (86, -4),
                (66, 0),
                (90, 0),
                (88, -1),
                (42, -5),
            ])
        ),
        (
            73,
            HashMap::from([
                (65, -2),
                (82, -4),
                (78, -4),
                (68, -4),
                (67, -2),
                (81, -4),
                (69, -4),
                (71, -5),
                (72, -4),
                (73, 6),
                (76, 2),
                (75, -4),
                (77, 2),
                (70, 0),
                (80, -3),
                (83, -3),
                (84, -1),
                (87, -3),
                (89, -1),
                (86, 4),
                (66, -4),
                (90, -4),
                (88, -1),
                (42, -5),
            ])
        ),
        (
            76,
            HashMap::from([
                (65, -2),
                (82, -3),
                (78, -4),
                (68, -5),
                (67, -2),
                (81, -3),
                (69, -4),
                (71, -5),
                (72, -3),
                (73, 2),
                (76, 6),
                (75, -3),
                (77, 3),
                (70, 1),
                (80, -4),
                (83, -3),
                (84, -2),
                (87, -3),
                (89, -1),
                (86, 1),
                (66, -4),
                (90, -3),
                (88, -1),
                (42, -5),
            ])
        ),
        (
            75,
            HashMap::from([
                (65, -1),
                (82, 3),
                (78, 0),
                (68, -1),
                (67, -4),
                (81, 2),
                (69, 1),
                (71, -2),
                (72, 0),
                (73, -4),
                (76, -3),
                (75, 6),
                (77, -2),
                (70, -4),
                (80, -1),
                (83, 0),
                (84, -1),
                (87, -4),
                (89, -2),
                (86, -3),
                (66, 0),
                (90, 1),
                (88, -1),
                (42, -5),
            ])
        ),
        (
            77,
            HashMap::from([
                (65, -1),
                (82, -2),
                (78, -3),
                (68, -4),
                (67, -2),
                (81, 0),
                (69, -3),
                (71, -3),
                (72, -2),
                (73, 2),
                (76, 3),
                (75, -2),
                (77, 8),
                (70, 0),
                (80, -3),
                (83, -2),
                (84, -1),
                (87, -2),
                (89, -1),
                (86, 1),
                (66, -3),
                (90, -2),
                (88, -1),
                (42, -5),
            ])
        ),
        (
            70,
            HashMap::from([
                (65, -3),
                (82, -3),
                (78, -4),
                (68, -5),
                (67, -3),
                (81, -4),
                (69, -4),
                (71, -4),
                (72, -1),
                (73, 0),
                (76, 1),
                (75, -4),
                (77, 0),
                (70, 9),
                (80, -5),
                (83, -3),
                (84, -3),
                (87, 2),
                (89, 4),
                (86, -1),
                (66, -5),
                (90, -4),
                (88, -2),
                (42, -5),
            ])
        ),
        (
            80,
            HashMap::from([
                (65, -1),
                (82, -3),
                (78, -2),
                (68, -2),
                (67, -3),
                (81, -1),
                (69, -1),
                (71, -3),
                (72, -3),
                (73, -3),
                (76, -4),
                (75, -1),
                (77, -3),
                (70, -5),
                (80, 10),
                (83, -1),
                (84, -1),
                (87, -5),
                (89, -4),
                (86, -3),
                (66, -2),
                (90, -1),
                (88, -2),
                (42, -5),
            ])
        ),
        (
            83,
            HashMap::from([
                (65, 2),
                (82, -1),
                (78, 1),
                (68, 0),
                (67, -1),
                (81, 0),
                (69, 0),
                (71, 0),
                (72, -1),
                (73, -3),
                (76, -3),
                (75, 0),
                (77, -2),
                (70, -3),
                (80, -1),
                (83, 5),
                (84, 2),
                (87, -4),
                (89, -2),
                (86, -2),
                (66, 0),
                (90, 0),
                (88, -1),
                (42, -5),
            ])
        ),
        (
            84,
            HashMap::from([
                (65, 0),
                (82, -1),
                (78, 0),
                (68, -1),
                (67, -1),
                (81, -1),
                (69, -1),
                (71, -2),
                (72, -2),
                (73, -1),
                (76, -2),
                (75, -1),
                (77, -1),
                (70, -3),
                (80, -1),
                (83, 2),
                (84, 6),
                (87, -3),
                (89, -2),
                (86, 0),
                (66, -1),
                (90, -1),
                (88, -1),
                (42, -5),
            ])
        ),
        (
            87,
            HashMap::from([
                (65, -4),
                (82, -3),
                (78, -5),
                (68, -5),
                (67, -4),
                (81, -2),
                (69, -3),
                (71, -3),
                (72, -3),
                (73, -3),
                (76, -3),
                (75, -4),
                (77, -2),
                (70, 2),
                (80, -5),
                (83, -4),
                (84, -3),
                (87, 15),
                (89, 3),
                (86, -4),
                (66, -5),
                (90, -3),
                (88, -3),
                (42, -5),
            ])
        ),
        (
            89,
            HashMap::from([
                (65, -2),
                (82, -2),
                (78, -2),
                (68, -3),
                (67, -3),
                (81, -1),
                (69, -2),
                (71, -4),
                (72, 2),
                (73, -1),
                (76, -1),
                (75, -2),
                (77, -1),
                (70, 4),
                (80, -4),
                (83, -2),
                (84, -2),
                (87, 3),
                (89, 9),
                (86, -2),
                (66, -3),
                (90, -2),
                (88, -1),
                (42, -5),
            ])
        ),
        (
            86,
            HashMap::from([
                (65, 0),
                (82, -3),
                (78, -4),
                (68, -4),
                (67, -1),
                (81, -3),
                (69, -3),
                (71, -4),
                (72, -4),
                (73, 4),
                (76, 1),
                (75, -3),
                (77, 1),
                (70, -1),
                (80, -3),
                (83, -2),
                (84, 0),
                (87, -4),
                (89, -2),
                (86, 5),
                (66, -4),
                (90, -3),
                (88, -1),
                (42, -5),
            ])
        ),
        (
            66,
            HashMap::from([
                (65, -2),
                (82, -1),
                (78, 4),
                (68, 5),
                (67, -4),
                (81, 0),
                (69, 1),
                (71, -1),
                (72, 0),
                (73, -4),
                (76, -4),
                (75, 0),
                (77, -3),
                (70, -5),
                (80, -2),
                (83, 0),
                (84, -1),
                (87, -5),
                (89, -3),
                (86, -4),
                (66, 5),
                (90, 2),
                (88, -1),
                (42, -5),
            ])
        ),
        (
            90,
            HashMap::from([
                (65, -1),
                (82, 0),
                (78, 0),
                (68, 1),
                (67, -4),
                (81, 4),
                (69, 5),
                (71, -3),
                (72, 0),
                (73, -4),
                (76, -3),
                (75, 1),
                (77, -2),
                (70, -4),
                (80, -1),
                (83, 0),
                (84, -1),
                (87, -3),
                (89, -2),
                (86, -3),
                (66, 2),
                (90, 5),
                (88, -1),
                (42, -5),
            ])
        ),
        (
            88,
            HashMap::from([
                (65, -1),
                (82, -1),
                (78, -1),
                (68, -2),
                (67, -2),
                (81, -1),
                (69, -1),
                (71, -2),
                (72, -1),
                (73, -1),
                (76, -1),
                (75, -1),
                (77, -1),
                (70, -2),
                (80, -2),
                (83, -1),
                (84, -1),
                (87, -3),
                (89, -1),
                (86, -1),
                (66, -1),
                (90, -1),
                (88, -1),
                (42, -5),
            ])
        ),
        (
            42,
            HashMap::from([
                (65, -5),
                (82, -5),
                (78, -5),
                (68, -5),
                (67, -5),
                (81, -5),
                (69, -5),
                (71, -5),
                (72, -5),
                (73, -5),
                (76, -5),
                (75, -5),
                (77, -5),
                (70, -5),
                (80, -5),
                (83, -5),
                (84, -5),
                (87, -5),
                (89, -5),
                (86, -5),
                (66, -5),
                (90, -5),
                (88, -5),
                (42, 1),
            ])
        ),
    ]);
}
