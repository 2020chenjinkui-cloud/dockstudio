// @vitest-environment jsdom
import { cleanup, fireEvent, render, screen } from '@testing-library/react';
import { afterEach, describe, expect, it } from 'vitest';
import FormatHint from './FormatHint';

afterEach(cleanup);

describe('FormatHint', () => {
  it('默认收起，只显示"格式与数据来源"一行', () => {
    render(<FormatHint kind="receptor" />);
    expect(screen.getByText('格式与数据来源')).toBeTruthy();
    expect(screen.queryByText(/受体文件/)).toBeNull();
  });

  it('点击展开内容，再点收起', () => {
    render(<FormatHint kind="receptor" />);
    const toggle = screen.getByRole('button');
    expect(toggle.getAttribute('aria-expanded')).toBe('false');
    fireEvent.click(toggle);
    expect(toggle.getAttribute('aria-expanded')).toBe('true');
    expect(screen.getByText(/受体文件/)).toBeTruthy();
    expect(screen.getByText(/rcsb\.org\/structure/)).toBeTruthy();
    fireEvent.click(toggle);
    expect(screen.queryByText(/受体文件/)).toBeNull();
  });

  it('library 变体包含 .smi 示例', () => {
    render(<FormatHint kind="library" />);
    fireEvent.click(screen.getByRole('button'));
    expect(screen.getByText(/CCO ethanol/)).toBeTruthy();
    expect(screen.getByText(/smiles 列/)).toBeTruthy();
  });

  it('ligand 变体包含盐拆分提醒', () => {
    render(<FormatHint kind="ligand" />);
    fireEvent.click(screen.getByRole('button'));
    expect(screen.getByText(/盐类/)).toBeTruthy();
    expect(screen.getByText(/PubChem/)).toBeTruthy();
  });
});
