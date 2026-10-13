# Chiffres complémentaires pour la révision

Les résultats de v1.0.0 sont conservés. Les analyses complémentaires utilisent le même jeu analytique et des erreurs standards clusterisées par paire pour les modèles appariés.

## Fibrose

Sur les mêmes 372 paires, le modèle complet sans fibrose donne +2.26 (IC95 % -0.38 à 4.90 ; p = 0.093) mm²/min, contre +2.36 (IC95 % -0.26 à 4.99 ; p = 0.078) avec fibrose. Dans les grandes lésions, le modèle année + taille + fibrose donne +5.29 (IC95 % 0.86 à 9.71 ; p = 0.019) mm²/min.

La restriction seule, avec le modèle année + taille, donne +3.00 (IC95 % 0.27 à 5.72 ; p = 0.031) mm²/min. Avec toutes les covariables et sans fibrose, les 391 paires initiales donnent +1.97 (IC95 % -0.58 à 4.51 ; p = 0.130). L'atténuation n'est donc pas expliquée par l'ajout de la fibrose. Elle est associée à l'ajustement supplémentaire ; la seule comparaison des modèles ne permet pas d'attribuer cette différence à la localisation ou à la restriction exclusivement.

## Taille continue

Il s'agit d'un test de Wald conjoint, avec covariance robuste clusterisée par paire, sur les quatre coefficients d'interaction entre ATRACT et les quatre fonctions de base du B-spline de taille (χ² = 12.842, 4 degrés de liberté ; p = 0.012073). La courbe est numériquement identique à table_s9 de v1.0.0. Les intervalles affichés sont ponctuels à 95 %, non simultanés.

## Volume des opérateurs

La classification est issue de la définition clinique fournie le 13 octobre 2026 : operator_01 et operator_02, environ 300 ESD/an et plus de 10 ans d'expérience ; les autres, tous formés, environ 50–100 ESD/an et environ 5 ans d'expérience. Cette classification fixe n'est pas une mesure annuelle de l'expérience au moment de chaque geste. Le modèle d'interaction ajuste sur les covariables préopératoires, y compris les effets fixes des opérateurs, sans fibrose.

lower-volume trained operators : n = 196, -2.08 (IC95 % -5.79 à 1.63 ; p = 0.271) mm²/min.

high-volume operators : n = 586, +3.40 (IC95 % 0.18 à 6.63 ; p = 0.039) mm²/min.

operator-volume interaction : n = 782, +5.49 (IC95 % 0.48 à 10.50 ; p = 0.032) mm²/min.

## Durée directement mesurée

overall : différence de durée ATRACT moins traction conventionnelle = -4.72 (IC95 % -9.97 à 0.53 ; p = 0.078) minutes. Une valeur négative correspond à une durée plus courte.

>=50 mm : différence de durée ATRACT moins traction conventionnelle = -9.65 (IC95 % -18.32 à -0.97 ; p = 0.029) minutes. Une valeur négative correspond à une durée plus courte.

Le modèle ajuste sur l'année et un B-spline à quatre degrés de liberté de la surface du spécimen. Les effets par taille proviennent d'un seul modèle d'interaction. Cette estimation de durée totale est distincte d'une conversion de vitesse en minutes.

## Calipers élargis

Caliper 0.40, overall : 436 paires, n = 872, +3.04 (IC95 % 0.39 à 5.68 ; p = 0.024) mm²/min ; SMD maximal 0.252.

Caliper 0.40, >=50 mm : 436 paires, n = 467, +5.97 (IC95 % 1.54 à 10.40 ; p = 0.008) mm²/min ; SMD maximal 0.252.

Caliper 0.40, size interaction : 436 paires, n = 872, +6.17 (IC95 % 1.00 à 11.34 ; p = 0.019) mm²/min ; SMD maximal 0.252.

Caliper 0.60, overall : 465 paires, n = 930, +3.30 (IC95 % 0.69 à 5.92 ; p = 0.013) mm²/min ; SMD maximal 0.381.

Caliper 0.60, >=50 mm : 465 paires, n = 503, +6.40 (IC95 % 2.07 à 10.73 ; p = 0.004) mm²/min ; SMD maximal 0.381.

Caliper 0.60, size interaction : 465 paires, n = 930, +6.55 (IC95 % 1.34 à 11.76 ; p = 0.014) mm²/min ; SMD maximal 0.381.

Les calipers 0,40 et 0,60 ne satisfont pas le critère de balance SMD <0,10. Ils restent des analyses exploratoires de relâchement, sans remplacer l'analyse principale ni permettre de conclure à une meilleure estimation causale.

## Sources

Toutes les estimations sont dans results/tables/estimates_summary.csv. Les effectifs n correspondent aux procédures de la strate ; n_model et n_pairs_model donnent la population utilisée pour estimer le modèle. Les termes d'interaction correspondent à une différence entre effets, non à un effet dans une strate supplémentaire. Aucun test n'est sélectionné selon sa significativité.
