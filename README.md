![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

# GenBank-to-BioDataset  
Pipeline de data engineering + bioinformatique pour produire des datasets ML-ready à partir de NCBI GenBank.

## 1. Titre & tagline
GenBank-to-BioDataset – Transformer des requêtes GenBank ciblées en jeux de données exploitables pour l’IA et la bioinformatique, avec un module de consensus expérimental pour les gènes barcodes (COI/COX1).

## 2. Description du projet
Les données publiques GenBank sont riches mais rarement prêtes pour l’entraînement de modèles ou les analyses reproductibles. GenBank-to-BioDataset propose une approche data engineering + bioinformatique : chaque requête est convertie à la demande en un pipeline complet (E-utilities uniquement), respectueux des bonnes pratiques NCBI et orienté vers les besoins IA/stage/medtech.

## 3. Fonctionnalités principales
- Téléchargement automatisé (esearch/efetch/esummary) puis export FASTA/CSV/JSONL.
- Nettoyage et filtrage qualité (A/C/G/T/N, longueur, fraction de N).
- Extraction de features ADN (longueur, %GC, comptages de bases).
- Rate limiting configuré et cache local hashé pour éviter les requêtes redondantes.
- CLI reproductible (`genbank2bio run`) avec configuration YAML.
- Module expérimental `genbank2bio consensus` centré sur COI/COX1 pour générer un consensus probable par clustering d’amplicons partiels.

## 4. Installation (Windows 10/11 – Python 3.10+)
```powershell
git clone https://github.com/you/genbank-to-biodataset.git
cd genbank-to-biodataset
py -3.10 -m venv .venv
.\.venv\Scripts\activate
pip install --upgrade pip
pip install -e .
```
Configurer `config/config.yaml` à partir de `config/example.yaml` (email, tool, rate_limit, etc.).

## 5. Démo rapide – Dataset GenBank
```powershell
genbank2bio run --query "Mammuthus primigenius[Organism] AND mitochondrion[Filter] AND complete genome" --limit 20
```
Sorties (par défaut `data/processed/`) :
- `sequences.fasta` : séquences nettoyées (FASTA).
- `features.csv` : features ADN prêtes pour un modèle (longueur, %GC, counts A/C/G/T/N).
- `metadata.jsonl` : métadonnées légères (organism, titre, dates, longueur).
Le résumé CLI indique le nombre d’IDs récupérés et les chemins exacts des fichiers.

## 6. Démo avancée – Consensus COI (expérimental)
**Use case** : COI/COX1 est un gène barcode standard pour différencier des espèces animales. Les séquences publiques sont souvent fragmentées ; le module consensus tente une reconstruction probable par clustering de fragments.

Commande :
```powershell
genbank2bio consensus --query "Aves[Organism] AND (COI[Gene] OR COX1[Title] OR \"cytochrome c oxidase subunit I\"[Title]) AND mitochondrion[Filter]" --limit 200
```
Sorties (`data/processed_consensus/` par défaut) :
- `clusters.jsonl` : contenu des clusters (accessions, tailles, statistiques).
- `consensus.fasta` : un consensus expérimental par cluster (header = taille cluster, longueur, %N, support moyen).
- `consensus_stats.csv` : métriques par cluster (n_seqs, longueur, %GC, fraction N, support).
- `filtered_features.csv` : features individuelles des séquences retenues après filtrage.
- `filtered_metadata.jsonl` : métadonnées associées aux séquences filtrées.

## 7. Interprétation scientifique & limites
- Les consensus générés ne constituent pas une vérité biologique ; ils représentent une hypothèse basée sur les fragments disponibles.
- Sensibilité forte au clustering, aux paramètres de filtrage et aux biais d’échantillonnage.
- Usage recommandé : exploratoire, pédagogique, prototypage IA. Toute utilisation clinique ou réglementaire requiert des validations supplémentaires.
- Validation indispensable : alignement manuel, comparaison à une référence, BLAST (prévu en V2).

## 8. Respect des règles NCBI
- Pas de mirroring massif : seules des requêtes ciblées via E-utilities sont exécutées.
- Rate limit configurable (0,35 s par défaut) pour respecter la politique NCBI.
- Cache local hashé pour limiter les re-téléchargements.
- L’utilisateur fournit son email/identité d’outil (champ `tool`) comme demandé par NCBI.

## 9. Roadmap (V2)
- Ajout optionnel de BLAST pour valider/annoter les consensus et enrichir les métadonnées.
- Amélioration des algorithmes de clustering (vsearch, heuristiques k-mer avancées).
- Support taxonomique plus fin (filtres hiérarchiques, métadonnées additionnelles).
- Nouvelles features ADN (motifs, entropie, signatures phylogénétiques).

## 10. Public visé / cas d’usage
- Candidats et encadrants de stage IA/bioinformatique souhaitant démontrer un workflow complet.
- Laboratoires ou équipes académiques voulant prototyper rapidement des pipelines GenBank-to-ML.
- Startups medtech/healthtech qui explorent des approches ML sur données génomiques publiques.
- Formations ou cours qui ont besoin d’un support reproductible pour illustrer data engineering + bioinfo.

## 11. Licence & citation
- Licence : MIT (voir fichier LICENCE si présent).
- Sources de données : NCBI GenBank (veiller à citer NCBI et les publications originales lors d’un usage scientifique ou public).
