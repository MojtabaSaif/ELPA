#include"Type.h"
#include<iostream>

vector<PACKET>::iterator NodeLabeles::Find(LABEL_ID id)
{
    if (size()==1)
		return (begin()->first == id) ? begin() : end();
	
	vector<PACKET>::iterator it = begin();
	for (; it != end(); it++)
		if (it->first == id)
			return it;
	
	return end();
}
