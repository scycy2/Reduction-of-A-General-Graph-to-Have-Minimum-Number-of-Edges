//
//  Rect.h
//  forceModel
//
//  Created by Ycl's Pro on 2021/6/12.
//

#ifndef Rect_h
#define Rect_h

class Rect {
private:
    int top, bottom, right, left;
public:
    Rect();
    Rect(int l, int t, int r, int b);
    ~Rect(){};
    void Assign(int l, int t, int r, int b);
};


#endif /* Rect_h */
