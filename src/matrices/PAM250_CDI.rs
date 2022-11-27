use super::Matrix;
use std::collections::HashMap;

lazy_static! {
    pub static ref MATRIX: Matrix = HashMap::from([
        (
            65,
            HashMap::from([
                (65, 6),
                (82, -5),
                (78, 1),
                (68, 1),
                (67, -7),
                (81, -1),
                (69, 1),
                (71, 4),
                (72, -4),
                (73, -2),
                (76, -6),
                (75, -4),
                (77, -4),
                (70, -12),
                (80, 4),
                (83, 4),
                (84, 4),
                (87, -19),
                (89, -12),
                (86, 1),
                (66, 1),
                (90, 0),
                (88, -1),
                (42, -26),
            ])
        ),
        (
            82,
            HashMap::from([
                (65, -5),
                (82, 20),
                (78, 0),
                (68, -4),
                (67, -12),
                (81, 4),
                (69, -4),
                (71, -8),
                (72, 5),
                (73, -7),
                (76, -10),
                (75, 11),
                (77, -1),
                (70, -15),
                (80, -1),
                (83, -1),
                (84, -3),
                (87, 7),
                (89, -14),
                (86, -8),
                (66, -2),
                (90, 0),
                (88, -3),
                (42, -26),
            ])
        ),
        (
            78,
            HashMap::from([
                (65, 1),
                (82, 0),
                (78, 7),
                (68, 7),
                (67, -12),
                (81, 3),
                (69, 5),
                (71, 1),
                (72, 5),
                (73, -6),
                (76, -10),
                (75, 3),
                (77, -6),
                (70, -12),
                (80, -2),
                (83, 2),
                (84, 1),
                (87, -13),
                (89, -7),
                (86, -6),
                (66, 7),
                (90, 4),
                (88, -1),
                (42, -26),
            ])
        ),
        (
            68,
            HashMap::from([
                (65, 1),
                (82, -4),
                (78, 7),
                (68, 13),
                (67, -17),
                (81, 5),
                (69, 11),
                (71, 2),
                (72, 2),
                (73, -8),
                (76, -13),
                (75, 0),
                (77, -9),
                (70, -19),
                (80, -3),
                (83, 1),
                (84, 0),
                (87, -22),
                (89, -14),
                (86, -7),
                (66, 10),
                (90, 9),
                (88, -2),
                (42, -26),
            ])
        ),
        (
            67,
            HashMap::from([
                (65, -7),
                (82, -12),
                (78, -12),
                (68, -17),
                (67, 40),
                (81, -18),
                (69, -18),
                (71, -11),
                (72, -11),
                (73, -8),
                (76, -20),
                (75, -18),
                (77, -17),
                (70, -14),
                (80, -9),
                (83, 0),
                (84, -7),
                (87, -26),
                (89, 1),
                (86, -6),
                (66, -15),
                (90, -18),
                (88, -10),
                (42, -26),
            ])
        ),
        (
            81,
            HashMap::from([
                (65, -1),
                (82, 4),
                (78, 3),
                (68, 5),
                (67, -18),
                (81, 13),
                (69, 8),
                (71, -4),
                (72, 10),
                (73, -7),
                (76, -6),
                (75, 2),
                (77, -3),
                (70, -15),
                (80, 1),
                (83, -2),
                (84, -3),
                (87, -16),
                (89, -13),
                (86, -6),
                (66, 4),
                (90, 11),
                (88, -2),
                (42, -26),
            ])
        ),
        (
            69,
            HashMap::from([
                (65, 1),
                (82, -4),
                (78, 5),
                (68, 11),
                (67, -18),
                (81, 8),
                (69, 13),
                (71, 1),
                (72, 2),
                (73, -7),
                (76, -11),
                (75, 0),
                (77, -7),
                (70, -18),
                (80, -2),
                (83, 0),
                (84, -1),
                (87, -23),
                (89, -14),
                (86, -6),
                (66, 9),
                (90, 11),
                (88, -2),
                (42, -26),
            ])
        ),
        (
            71,
            HashMap::from([
                (65, 4),
                (82, -8),
                (78, 1),
                (68, 2),
                (67, -11),
                (81, -4),
                (69, 1),
                (71, 16),
                (72, -7),
                (73, -8),
                (76, -13),
                (75, -6),
                (77, -9),
                (70, -16),
                (80, -2),
                (83, 4),
                (84, 0),
                (87, -23),
                (89, -17),
                (86, -4),
                (66, 2),
                (90, -1),
                (88, -3),
                (42, -26),
            ])
        ),
        (
            72,
            HashMap::from([
                (65, -4),
                (82, 5),
                (78, 5),
                (68, 2),
                (67, -11),
                (81, 10),
                (69, 2),
                (71, -7),
                (72, 22),
                (73, -8),
                (76, -7),
                (75, 0),
                (77, -7),
                (70, -6),
                (80, -1),
                (83, -3),
                (84, -4),
                (87, -9),
                (89, 0),
                (86, -7),
                (66, 4),
                (90, 6),
                (88, -2),
                (42, -26),
            ])
        ),
        (
            73,
            HashMap::from([
                (65, -2),
                (82, -7),
                (78, -6),
                (68, -8),
                (67, -8),
                (81, -7),
                (69, -7),
                (71, -8),
                (72, -8),
                (73, 15),
                (76, 8),
                (75, -6),
                (77, 7),
                (70, 3),
                (80, -7),
                (83, -5),
                (84, 0),
                (87, -17),
                (89, -3),
                (86, 12),
                (66, -7),
                (90, -7),
                (88, -2),
                (42, -26),
            ])
        ),
        (
            76,
            HashMap::from([
                (65, -6),
                (82, -10),
                (78, -10),
                (68, -13),
                (67, -20),
                (81, -6),
                (69, -11),
                (71, -13),
                (72, -7),
                (73, 8),
                (76, 20),
                (75, -9),
                (77, 12),
                (70, 6),
                (80, -8),
                (83, -9),
                (84, -6),
                (87, -6),
                (89, -3),
                (86, 6),
                (66, -11),
                (90, -9),
                (88, -4),
                (42, -26),
            ])
        ),
        (
            75,
            HashMap::from([
                (65, -4),
                (82, 11),
                (78, 3),
                (68, 0),
                (67, -18),
                (81, 2),
                (69, 0),
                (71, -6),
                (72, 0),
                (73, -6),
                (76, -9),
                (75, 15),
                (77, 1),
                (70, -17),
                (80, -4),
                (83, -1),
                (84, 0),
                (87, -12),
                (89, -15),
                (86, -8),
                (66, 2),
                (90, 1),
                (88, -3),
                (42, -26),
            ])
        ),
        (
            77,
            HashMap::from([
                (65, -4),
                (82, -1),
                (78, -6),
                (68, -9),
                (67, -17),
                (81, -3),
                (69, -7),
                (71, -9),
                (72, -7),
                (73, 7),
                (76, 12),
                (75, 1),
                (77, 21),
                (70, 1),
                (80, -7),
                (83, -5),
                (84, -2),
                (87, -14),
                (89, -8),
                (86, 6),
                (66, -7),
                (90, -5),
                (88, -2),
                (42, -26),
            ])
        ),
        (
            70,
            HashMap::from([
                (65, -12),
                (82, -15),
                (78, -12),
                (68, -19),
                (67, -14),
                (81, -15),
                (69, -18),
                (71, -16),
                (72, -6),
                (73, 3),
                (76, 6),
                (75, -17),
                (77, 1),
                (70, 30),
                (80, -15),
                (83, -11),
                (84, -10),
                (87, 1),
                (89, 23),
                (86, -4),
                (66, -15),
                (90, -17),
                (88, -8),
                (42, -26),
            ])
        ),
        (
            80,
            HashMap::from([
                (65, 4),
                (82, -1),
                (78, -2),
                (68, -3),
                (67, -9),
                (81, 1),
                (69, -2),
                (71, -2),
                (72, -1),
                (73, -7),
                (76, -8),
                (75, -4),
                (77, -7),
                (70, -15),
                (80, 19),
                (83, 3),
                (84, 1),
                (87, -19),
                (89, -16),
                (86, -4),
                (66, -2),
                (90, -1),
                (88, -2),
                (42, -26),
            ])
        ),
        (
            83,
            HashMap::from([
                (65, 4),
                (82, -1),
                (78, 2),
                (68, 1),
                (67, 0),
                (81, -2),
                (69, 0),
                (71, 4),
                (72, -3),
                (73, -5),
                (76, -9),
                (75, -1),
                (77, -5),
                (70, -11),
                (80, 3),
                (83, 5),
                (84, 4),
                (87, -8),
                (89, -9),
                (86, -3),
                (66, 2),
                (90, -1),
                (88, -1),
                (42, -26),
            ])
        ),
        (
            84,
            HashMap::from([
                (65, 4),
                (82, -3),
                (78, 1),
                (68, 0),
                (67, -7),
                (81, -3),
                (69, -1),
                (71, 0),
                (72, -4),
                (73, 0),
                (76, -6),
                (75, 0),
                (77, -2),
                (70, -10),
                (80, 1),
                (83, 4),
                (84, 9),
                (87, -17),
                (89, -9),
                (86, 1),
                (66, 0),
                (90, -2),
                (88, -1),
                (42, -26),
            ])
        ),
        (
            87,
            HashMap::from([
                (65, -19),
                (82, 7),
                (78, -13),
                (68, -22),
                (67, -26),
                (81, -16),
                (69, -23),
                (71, -23),
                (72, -9),
                (73, -17),
                (76, -6),
                (75, -12),
                (77, -14),
                (70, 1),
                (80, -19),
                (83, -8),
                (84, -17),
                (87, 57),
                (89, 0),
                (86, -21),
                (66, -18),
                (90, -19),
                (88, -13),
                (42, -26),
            ])
        ),
        (
            89,
            HashMap::from([
                (65, -12),
                (82, -14),
                (78, -7),
                (68, -14),
                (67, 1),
                (81, -13),
                (69, -14),
                (71, -17),
                (72, 0),
                (73, -3),
                (76, -3),
                (75, -15),
                (77, -8),
                (70, 23),
                (80, -16),
                (83, -9),
                (84, -9),
                (87, 0),
                (89, 34),
                (86, -8),
                (66, -10),
                (90, -14),
                (88, -8),
                (42, -26),
            ])
        ),
        (
            86,
            HashMap::from([
                (65, 1),
                (82, -8),
                (78, -6),
                (68, -7),
                (67, -6),
                (81, -6),
                (69, -6),
                (71, -4),
                (72, -7),
                (73, 12),
                (76, 6),
                (75, -8),
                (77, 6),
                (70, -4),
                (80, -4),
                (83, -3),
                (84, 1),
                (87, -21),
                (89, -8),
                (86, 14),
                (66, -6),
                (90, -6),
                (88, -2),
                (42, -26),
            ])
        ),
        (
            66,
            HashMap::from([
                (65, 1),
                (82, -2),
                (78, 7),
                (68, 10),
                (67, -15),
                (81, 4),
                (69, 9),
                (71, 2),
                (72, 4),
                (73, -7),
                (76, -11),
                (75, 2),
                (77, -7),
                (70, -15),
                (80, -2),
                (83, 2),
                (84, 0),
                (87, -18),
                (89, -10),
                (86, -6),
                (66, 9),
                (90, 7),
                (88, -2),
                (42, -26),
            ])
        ),
        (
            90,
            HashMap::from([
                (65, 0),
                (82, 0),
                (78, 4),
                (68, 9),
                (67, -18),
                (81, 11),
                (69, 11),
                (71, -1),
                (72, 6),
                (73, -7),
                (76, -9),
                (75, 1),
                (77, -5),
                (70, -17),
                (80, -1),
                (83, -1),
                (84, -2),
                (87, -19),
                (89, -14),
                (86, -6),
                (66, 7),
                (90, 11),
                (88, -2),
                (42, -26),
            ])
        ),
        (
            88,
            HashMap::from([
                (65, -1),
                (82, -3),
                (78, -1),
                (68, -2),
                (67, -10),
                (81, -2),
                (69, -2),
                (71, -3),
                (72, -2),
                (73, -2),
                (76, -4),
                (75, -3),
                (77, -2),
                (70, -8),
                (80, -2),
                (83, -1),
                (84, -1),
                (87, -13),
                (89, -8),
                (86, -2),
                (66, -2),
                (90, -2),
                (88, -3),
                (42, -26),
            ])
        ),
        (
            42,
            HashMap::from([
                (65, -26),
                (82, -26),
                (78, -26),
                (68, -26),
                (67, -26),
                (81, -26),
                (69, -26),
                (71, -26),
                (72, -26),
                (73, -26),
                (76, -26),
                (75, -26),
                (77, -26),
                (70, -26),
                (80, -26),
                (83, -26),
                (84, -26),
                (87, -26),
                (89, -26),
                (86, -26),
                (66, -26),
                (90, -26),
                (88, -26),
                (42, 1),
            ])
        ),
    ]);
}