#!/usr/bin/env python
# coding: utf-8

# In[1]:


print("Hello World")


# In[3]:


r = 1
E_r = {
    (0, 0): VectorSpace(QQ, 1),  # 1-dimensional vector space
    (1, 0): VectorSpace(QQ, 2),  # 2-dimensional vector space
    (0, 1): VectorSpace(QQ, 1),  # 1-dimensional vector space
    (1, 1): VectorSpace(QQ, 1),  # 1-dimensional vector space
}

# Example: Define differentials d_r
differentials = {
    ((0, 0), (1, 0)): matrix(QQ, [[1, 0]]),  # d_r: E_r^{0,0} -> E_r^{1,0}
    ((0, 1), (1, 1)): matrix(QQ, [[1]]),     # d_r: E_r^{0,1} -> E_r^{1,1}
}

def compute_next_page(E_r, differentials):
    E_next = {}
    for (p, q) in E_r.keys():
        # Compute the kernel of d_r^{p,q}
        if ((p, q), (p + r, q - r + 1)) in differentials:
            d = differentials[((p, q), (p + r, q - r + 1))]
            ker = d.kernel()
        else:
            ker = E_r[(p, q)]

        # Compute the image of d_r^{p-r, q+r-1}
        if ((p - r, q + r - 1), (p, q)) in differentials:
            d_prev = differentials[((p - r, q + r - 1), (p, q))]
            im = d_prev.image()
        else:
            im = VectorSpace(QQ, 0)  # Zero space

        # Compute the homology
        E_next[(p, q)] = ker / im

    return E_next

# Example: Compute the next page
E_next = compute_next_page(E_r, differentials)

# Print the next page
print("Next page of the spectral sequence:")
for (p, q) in E_next.keys():
    print(f"E_{r+1}^{p,q} = {E_next[(p, q)]}")


# In[ ]:




