/* stopping_srim.h -- pouvoirs d'arret electroniques SRIM-2013 / SRIM-2008 / TRIM-98 pour iradina
 *   stopping_model = 3 : SRIM 2012 (SCOEF03)
 *   stopping_model = 4 : SRIM 2008 (SCOEF03, meme fichier)
 *   stopping_model = 5 : TRIM 1998 (SCOEF.95A, SCOEF.95B, SCOEFGAS.95)
 * Les fonctions de calcul sont dans stop2012.c et stop1998.c ; ce module ne fait que
 * construire des tables par materiau et les interroger (meme convention que dpass). */
#ifndef STOPPING_SRIM_H
#define STOPPING_SRIM_H

/* Options par materiau, a appeler AVANT PrepareStoppingTablesSRIM (defaut : solide, compound=1).
   gas      : 1 = cible gazeuse (SRIM : coefficients gaz ; TRIM-98 : SCOEFGAS.95, necessaire alors)
   compound : facteur de correction de compose (modeles 3 et 4 ; <=0 equivaut a 1).
   Le modele 3 (SRIM-2013) a en plus la correction de He dans H ou C (gaz, compound != 1). */
void SetSRIMMaterialOptions(int material, int gas, double compound);

/* Lit les fichiers de donnees (dans DirectoryData) et construit, pour chaque materiau et
   chaque projectile possible (existing_elements), la table dE/dx [eV/nm].
   model = 3, 4 ou 5. Retourne 0 si OK, -1 memoire, -2 lecture des fichiers. */
int PrepareStoppingTablesSRIM(int model);

/* dE/dx electronique [eV/nm] de l'ion (Z, M en u, energie ionE en eV) dans le materiau. */
float ElectronicStoppingSRIM(int ionZ, float ionM, float ionE, int material);

#endif
