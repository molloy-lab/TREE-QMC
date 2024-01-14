#ifndef TOMS_H
#define TOMS_H

namespace toms743 {

double bisect ( double xx, int nb, int &ner, int l );
double crude ( double xx, int nb );
int nbits_compute ( );
void timestamp ( );
double wapr ( double x, int nb, int &nerror, int l );

}

#endif