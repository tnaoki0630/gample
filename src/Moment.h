#import <Foundation/Foundation.h>

@interface Moment : NSObject

struct Values {
    double n;
    double ux;
    double uy;
    double uz;
    double Pxx;
    double Pxy;
    double Pxz;
    double Pyy;
    double Pyz;
    double Pzz;
};

- (instancetype)initialize;

@end